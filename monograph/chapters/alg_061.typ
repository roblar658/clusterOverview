#import "../template/ieee.typ": *

#alg_entry(
  id: 61,
  name: "Centroid Linkage (UPGMC)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1958",
  paradigm: "Agglomerative Geometric Centroid Linkage",
  complexity: $cal(O)(n^2 log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $X = \{bold(x)_1, dots, bold(x)_n\} subset bb(R)^p$ be a finite point cloud in $p$-dimensional Euclidean space, equipped with the initial squared Euclidean dissimilarity matrix $D^((0))_(i j) = norm(bold(x)_i - bold(x)_j)_2^2$. Agglomerative hierarchical clustering produces a nested hierarchy of partitions $cal(C)^((0)), dots, cal(C)^((n-1))$ by greedily merging the closest pair of active clusters at each step $t = 1, dots, n - 1$.

Introduced by Sokal and Michener @sokal1958statistical as the *Unweighted Pair Group Method with Centroid* (*UPGMC*) and formulated combinatorially by Lance and Williams @lance1967general, Centroid Linkage represents each active cluster $A subset.eq X$ by its unweighted Euclidean barycenter (arithmetic mean) $bold(mu)_A = |A|^(-1) sum_(bold(x) in A) bold(x)$, and defines the dissimilarity between two disjoint clusters as the Euclidean distance between their centroids. Because each object inside $A$ contributes uniform weight $1/|A|$ to $bold(mu)_A$, the method is termed *unweighted* at the observation level, meaning that merging a large cluster $A$ with a small cluster $B$ shifts the combined centroid proportionally to $|B| / (|A| + |B|)$.

#definition(title: "UPGMC Centroid Linkage and Squared-Distance Recurrence")[
  For any non-empty cluster $A subset.eq X$, let $bold(mu)_A = |A|^(-1) sum_(bold(x) in A) bold(x)$ be its centroid. Define the Euclidean and squared-Euclidean centroid dissimilarities between disjoint clusters $A, B subset.eq X$ by
  $ d_C (A, B) = norm(bold(mu)_A - bold(mu)_B)_2, quad D_C (A, B) = norm(bold(mu)_A - bold(mu)_B)_2^2. $
  Starting from $cal(C)^((0)) = \{\{bold(x)_1\}, dots, \{bold(x)_n\}\}$, at step $t in \{1, dots, n - 1\}$ UPGMC selects a pair
  $ (A_t, B_t) in arg min_(U != V in cal(C)^((t-1))) d_C (U, V) = arg min_(U != V in cal(C)^((t-1))) D_C (U, V), $
  records merge height $h_t = d_C (A_t, B_t)$ (and squared height $H_t = h_t^2$), forms $A_t union B_t$ with centroid
  $ bold(mu)_(A_t union B_t) = (|A_t| bold(mu)_(A_t) + |B_t| bold(mu)_(B_t)) / (|A_t| + |B_t|), $
  and updates squared dissimilarities to every remaining cluster $C$ via the Lance--Williams recurrence ($alpha_A = |A|/(|A|+|B|)$, $alpha_B = |B|/(|A|+|B|)$, $beta = -|A||B|/(|A|+|B|)^2$, $gamma = 0$):
  $ D_C (A union B, C) &= (|A|) / (|A| + |B|) D_C (A, C) + (|B|) / (|A| + |B|) D_C (B, C) \ &quad - (|A| |B|) / ((|A| + |B|)^2) D_C (A, B). $
]

Unlike average or complete linkage, the negative coefficient $beta < 0$ breaks Bruynooghe reducibility, invalidating the $cal(O)(n^2)$ Nearest-Neighbor Chain algorithm @mullner2011modern. Deterministic exact execution therefore maintains per-cluster indexed min-heaps to achieve $cal(O)(n^2 log n)$ worst-case time.

#algbox(
  title: "Exact Centroid Linkage (UPGMC) via Indexed Heaps",
  inputs: [Points $X = \{bold(x)_1, dots, bold(x)_n\} subset bb(R)^p$ (or squared Euclidean matrix $D^((0))$).],
  outputs: [Merge list $\{(A_t, B_t, h_t)\}_(t=1)^(n-1)$ and cluster centroids.],
)[
  1. Compute $D_C (\{bold(x)_i\}, \{bold(x)_j\}) <- norm(bold(x)_i - bold(x)_j)_2^2$ for all $1 <= i < j <= n$; set $cal(C) <- \{\{bold(x)_1\}, dots, \{bold(x)_n\}\}$ and sizes $n_i <- 1$.
  2. For each active cluster $U in cal(C)$, build an indexed min-heap $Q_U$ storing $\{(D_C (U, V), V) : V in cal(C) without \{U\}\}$, and insert $(min(Q_U), U)$ into a global min-heap $G$.
  3. For $t = 1, dots, n - 1$:
     a. Extract the global minimizing pair $(A, B)$ with squared distance $H_t = D_C (A, B)$ from $G$; set $h_t <- sqrt(H_t)$ and record $(A, B, h_t)$.
     b. Create $W <- A union B$ with $n_W <- n_A + n_B$ and remove $A, B$ from $cal(C)$ and $G$.
     c. For each remaining $C in cal(C)$:
        - Delete $A$ and $B$ from $Q_C$.
        - Compute $D_C (W, C) <- (n_A / n_W) D_C (A, C) + (n_B / n_W) D_C (B, C) - (n_A n_B / n_W^2) H_t$.
        - Insert $(D_C (W, C), W)$ into $Q_C$, update $C$'s key in $G$, and collect $(D_C (W, C), C)$ for $Q_W$.
     d. Build heap $Q_W$ in $cal(O)(|cal(C)|)$ time, insert $(min(Q_W), W)$ into $G$, and add $W$ to $cal(C)$.
  4. Return the merge sequence $\{(A_t, B_t, h_t)\}_(t=1)^(n-1)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Hilbert-Space Stewart Identity and Combinatorial Exactness")

A foundational algebraic fact makes UPGMC executable directly from squared pairwise distances without ever storing or manipulating ambient coordinates in $bb(R)^p$.

#lemma(title: "Hilbert-Space Stewart Identity and Combinatoriality", origin: "Original Proof")[
  Let $A, B, C subset bb(R)^p$ be finite non-empty sets with $A inter B = emptyset$, and let $n_A = |A|$, $n_B = |B|$, $n_(A B) = n_A + n_B$. Then the squared Euclidean centroid distance satisfies
  $ D_C (A union B, C) = n_A / n_(A B) D_C (A, C) + n_B / n_(A B) D_C (B, C) - (n_A n_B) / n_(A B)^2 D_C (A, B). $
  Consequently, UPGMC is combinatorial on squared Euclidean distance matrices: every merge height and topology is uniquely determined by $D^((0))_(i j) = norm(bold(x)_i - bold(x)_j)_2^2$ alone, independently of $p$.
]

#proof[
  Set $bold(u) = bold(mu)_A - bold(mu)_C$, $bold(v) = bold(mu)_B - bold(mu)_C$, and $alpha = n_A / n_(A B) in (0, 1)$, so $1 - alpha = n_B / n_(A B)$. Then $bold(mu)_(A union B) - bold(mu)_C = alpha bold(u) + (1 - alpha) bold(v)$ and $bold(mu)_A - bold(mu)_B = bold(u) - bold(v)$. Expanding the inner product gives
  $ norm(alpha bold(u) + (1 - alpha) bold(v))_2^2 &= alpha^2 norm(bold(u))_2^2 + (1 - alpha)^2 norm(bold(v))_2^2 \ &quad + 2 alpha (1 - alpha) chevron.l bold(u), bold(v) chevron.r. $
  Substituting the polarization identity $2 chevron.l bold(u), bold(v) chevron.r = norm(bold(u))_2^2 + norm(bold(v))_2^2 - norm(bold(u) - bold(v))_2^2$ yields
  $ norm(alpha bold(u) + (1 - alpha) bold(v))_2^2 &= (alpha^2 + alpha (1 - alpha)) norm(bold(u))_2^2 \ &quad + ((1 - alpha)^2 + alpha (1 - alpha)) norm(bold(v))_2^2 \ &quad - alpha (1 - alpha) norm(bold(u) - bold(v))_2^2 \ &= alpha norm(bold(u))_2^2 + (1 - alpha) norm(bold(v))_2^2 \ &quad - alpha (1 - alpha) norm(bold(u) - bold(v))_2^2. $
  Replacing $alpha = n_A / n_(A B)$, $1 - alpha = n_B / n_(A B)$, $norm(bold(u))_2^2 = D_C (A, C)$, $norm(bold(v))_2^2 = D_C (B, C)$, and $norm(bold(u) - bold(v))_2^2 = D_C (A, B)$ yields the claimed formula. Since $D^((0))$ initializes $D_C$ for all singletons and cardinalities update via $n_(A B) = n_A + n_B$, induction on $t$ shows that coordinates in $bb(R)^p$ are never needed after computing $D^((0))$.
]

#subsec("B", "Dendrogram Inversions, Sharp Reversal Bound, and 1D Monotonicity")

The defining geometric pathology of UPGMC is the occurrence of *dendrogram inversions* (also called reversals), where $h_(t+1) < h_t$. We establish the exact dimensional threshold, the sharp single-step reversal lower bound, and a multi-step regular-simplex inversion chain.

#theorem(title: "Sharp Single-Step Inversion Bound, 1D Monotonicity, and Simplex Reversal Chain", origin: "Original Proof")[
  For any finite point set $X subset bb(R)^p$ clustered by UPGMC with merge heights $h_1, dots, h_(n-1)$:
  1. *Sharp single-step lower bound*: At every step $t in \{1, dots, n - 2\}$, merging $(A_t, B_t)$ yields
     $ h_(t+1) >= sqrt(1 - (|A_t| |B_t|) / (|A_t| + |B_t|)^2) h_t >= sqrt(3) / 2 h_t. $
  2. *Dimension threshold*: In $p = 1$ ($X subset bb(R)$), UPGMC is strictly inversion-free ($h_(t+1) >= h_t$ for all $t$). For every $p >= 2$, strict inversions occur and approach $h_(t+1) / h_t -> sqrt(3)/2$ arbitrarily closely; concretely, on $X = \{(-1, 0), (1, 0), (0, 1.8)\} subset bb(R)^2$, $h_1 = 2.00$ and $h_2 = 1.80 < h_1$.
  3. *Multi-step inversion chain*: On the $n$ vertices of a regular $(n - 1)$-simplex of unit edge length in $bb(R)^n$ (with infinitesimal tie-breaking), UPGMC undergoes $n - 2$ consecutive strict inversions with exact heights
     $ h_m = sqrt((m + 1) / (2 m)) arrow.b 1 / sqrt(2) quad (m = 1, dots, n - 1). $
]

#proof[
  For (1), at step $t$, $(A_t, B_t)$ minimizes $D_C$ over $cal(C)^((t-1))$, so $D_C (U, V) >= H_t = h_t^2$ for every untouched pair $(U, V)$, and $D_C (A_t, C) >= H_t$, $D_C (B_t, C) >= H_t$ for every remaining cluster $C$. Putting $alpha = |A_t| / (|A_t| + |B_t|) in (0, 1)$ in Lemma 1 gives
  $ D_C (A_t union B_t, C) &>= alpha H_t + (1 - alpha) H_t - alpha (1 - alpha) H_t \ &= (1 - alpha (1 - alpha)) H_t >= 3 / 4 H_t, $
  since $alpha (1 - alpha) <= 1/4$ with equality iff $|A_t| = |B_t|$. Taking the minimum over $cal(C)^((t))$ and square roots yields $h_(t+1) >= sqrt(1 - alpha (1 - alpha)) h_t >= (sqrt(3)/2) h_t$.

  For (2), let $p = 1$ and write $mu_(A_t) <= mu_(B_t)$ with $h_t = mu_(B_t) - mu_(A_t)$. For any remaining cluster $C$, if $mu_(A_t) < mu_C < mu_(B_t)$, then $(mu_C - mu_(A_t)) + (mu_(B_t) - mu_C) = h_t$, forcing $min(|mu_C - mu_(A_t)|, |mu_C - mu_(B_t)|) <= h_t / 2 < h_t$, contradicting minimality of $(A_t, B_t)$. Thus either $mu_C <= mu_(A_t) - h_t$ or $mu_C >= mu_(B_t) + h_t$ (or $h_t = 0$). Since $mu_(A_t union B_t) in [mu_(A_t), mu_(B_t)]$, in either case $|mu_C - mu_(A_t union B_t)| >= h_t$, proving $h_(t+1) >= h_t$. For $p >= 2$, take $bold(x)_1 = (-1, 0)$, $bold(x)_2 = (1, 0)$, $bold(x)_3 = (0, 1.8)$ (embedded in $bb(R)^p$): $d(bold(x)_1, bold(x)_2) = 2.00$, while $d(bold(x)_1, bold(x)_3) = d(bold(x)_2, bold(x)_3) = sqrt(1 + 3.24) = sqrt(4.24) approx 2.059 > 2.00$. Step 1 merges $A = \{bold(x)_1, bold(x)_2\}$ at $h_1 = 2.00$ with $bold(mu)_A = (0, 0)$, and Step 2 merges $A union \{bold(x)_3\}$ at $h_2 = 1.80 < h_1$. Replacing $1.8$ by $sqrt(3) + epsilon$ ($epsilon arrow.b 0$) drives $h_2 / h_1 -> sqrt(3)/2$.

  For (3), place $bold(x)_i = (1 / sqrt(2)) bold(e)_i in bb(R)^n$ ($i = 1, dots, n$), so $norm(bold(x)_i - bold(x)_j)_2 = 1$ for all $i != j$. With an infinitesimal perturbation breaking ties in favor of the growing cluster $S_m = \{bold(x)_1, dots, bold(x)_m\}$, its centroid is $bold(mu)_(S_m) = (m sqrt(2))^(-1) sum_(i=1)^m bold(e)_i$, orthogonal to every unmerged vertex $bold(x)_j$ ($j > m$). Hence
  $ D_C (S_m, \{bold(x)_j\}) = norm(bold(mu)_(S_m))_2^2 + norm(bold(x)_j)_2^2 = 1 / (2 m) + 1 / 2 = (m + 1) / (2 m), $
  which is strictly smaller than $D_C (\{bold(x)_j\}, \{bold(x)_k\}) = 1$ for all $m >= 2$ and strictly decreases as $m$ grows from $1$ to $n - 1$, converging to $1/2$ ($h_m arrow.b 1/sqrt(2)$).
]

#subsec("C", "Exact Variance Decomposition and Connection to Ward and UPGMA")

UPGMC admits two exact algebraic identities linking it directly to squared-Euclidean UPGMA (Alg. 59) and Ward's minimum-variance criterion (Alg. 58 @ward1963hierarchical).

#theorem(title: "Centroid-UPGMA Variance Gap and Ward Scaling Identity", origin: "Original Derivation")[
  For any finite cluster $A subset bb(R)^p$, let $"SSE"(A) = sum_(bold(x) in A) norm(bold(x) - bold(mu)_A)_2^2$ be its within-cluster sum of squares and $sigma_A^2 = |A|^(-1) "SSE"(A)$ its trace covariance. Let $D_U^((2)) (A, B) = (|A| |B|)^(-1) sum_(bold(x) in A, bold(y) in B) norm(bold(x) - bold(y))_2^2$ be squared-Euclidean UPGMA dissimilarity. Then for any disjoint clusters $A, B subset bb(R)^p$:
  1. *Centroid--UPGMA variance gap*:
     $ D_C (A, B) = D_U^((2)) (A, B) - sigma_A^2 - sigma_B^2. $
  2. *Ward--Centroid harmonic scaling*:
     $ Delta "SSE"(A, B) := "SSE"(A union B) - "SSE"(A) - "SSE"(B) = (|A| |B|) / (|A| + |B|) D_C (A, B). $
]

#proof[
  For (1), decompose $bold(x) - bold(y) = (bold(x) - bold(mu)_A) + (bold(mu)_A - bold(mu)_B) - (bold(y) - bold(mu)_B)$. Squaring and averaging uniformly over $(bold(x), bold(y)) in A times B$, all three cross inner products vanish identically because $sum_(bold(x) in A) (bold(x) - bold(mu)_A) = bold(0)$ and $sum_(bold(y) in B) (bold(y) - bold(mu)_B) = bold(0)$. Thus
  $ D_U^((2)) (A, B) &= 1 / |A| sum_(bold(x) in A) norm(bold(x) - bold(mu)_A)_2^2 + norm(bold(mu)_A - bold(mu)_B)_2^2 \ &quad + 1 / |B| sum_(bold(y) in B) norm(bold(y) - bold(mu)_B)_2^2 \ &= sigma_A^2 + D_C (A, B) + sigma_B^2. $
  For (2), by the parallel-axis theorem,
  $ "SSE"(A union B) &= "SSE"(A) + |A| norm(bold(mu)_A - bold(mu)_(A union B))_2^2 \ &quad + "SSE"(B) + |B| norm(bold(mu)_B - bold(mu)_(A union B))_2^2. $
  Since $bold(mu)_A - bold(mu)_(A union B) = (|B| / (|A| + |B|)) (bold(mu)_A - bold(mu)_B)$ and $bold(mu)_B - bold(mu)_(A union B) = -(|A| / (|A| + |B|)) (bold(mu)_A - bold(mu)_B)$, substituting yields
  $ Delta "SSE"(A, B) = (|A| |B|^2 + |B| |A|^2) / ((|A| + |B|)^2) D_C (A, B) = (|A| |B|) / (|A| + |B|) D_C (A, B). $
]

Identity (1) explains why UPGMC merges diffuse high-variance clusters earlier than squared-UPGMA (subtracting $sigma_A^2 + sigma_B^2$), while identity (2) shows UPGMC drops Ward's harmonic cardinality factor $|A||B|/(|A|+|B|)$, removing Ward's strong preference for equal-sized clusters.

#subsec("D", "Failure of the Nearest-Neighbor Chain Algorithm and Exact Complexity")

Because UPGMC violates Bruynooghe reducibility, eager reciprocal nearest-neighbor merging can return a provably incorrect dendrogram.

#proposition(title: "Constructive NN-Chain Failure and Deterministic Heap Complexity", origin: "Original Counterexample")[
  1. *NN-chain invalidity*: There exists a 4-point planar configuration $X subset bb(R)^2$ on which the Nearest-Neighbor Chain algorithm produces a strictly incorrect UPGMC dendrogram and incorrect 2-cluster partition @mullner2011modern.
  2. *Exact deterministic complexity*: Algorithm 61 (indexed per-cluster min-heaps) computes the exact UPGMC hierarchy in $cal(O)(n^2 log n)$ worst-case time and $cal(O)(n^2)$ space (plus $cal(O)(n^2 p)$ initial distance evaluation).
]

#proof[
  For (1), let $bold(x)_1 = (-1, 0)$, $bold(x)_2 = (1, 0)$, $bold(x)_3 = (0, 3.65)$, and $bold(x)_4 = (0, 1.80)$ in $bb(R)^2$. The pairwise Euclidean distances are $d(bold(x)_3, bold(x)_4) = 1.85$, $d(bold(x)_1, bold(x)_2) = 2.00$, $d(bold(x)_1, bold(x)_4) = d(bold(x)_2, bold(x)_4) = sqrt(4.24) approx 2.059$, and $d(bold(x)_1, bold(x)_3) = d(bold(x)_2, bold(x)_3) = sqrt(14.3225) approx 3.785$.
  - *Exact UPGMC*: Step 1 merges the global minimum $C_(3 4) = \{bold(x)_3, bold(x)_4\}$ at $h_1 = 1.85$ with centroid $(0, 2.725)$. Since $d(bold(x)_1, C_(3 4)) = sqrt(1 + 2.725^2) approx 2.903 > 2.00$, Step 2 merges $C_(1 2) = \{bold(x)_1, bold(x)_2\}$ at $h_2 = 2.00$, yielding the exact 2-cluster partition $\{\{bold(x)_1, bold(x)_2\}, \{bold(x)_3, bold(x)_4\}\}$.
  - *NN-chain initiated at $bold(x)_1$*: Since $2.00 < 2.059 < 3.785$, $"NN"(bold(x)_1) = bold(x)_2$ and $"NN"(bold(x)_2) = bold(x)_1$. The NN-chain eagerly merges the reciprocal pair $C_(1 2) = \{bold(x)_1, bold(x)_2\}$ with centroid $(0, 0)$. Immediately afterward, $d_C (C_(1 2), bold(x)_4) = 1.80 < 1.85 = d(bold(x)_3, bold(x)_4)$, causing the NN-chain to merge $C_(1 2) union \{bold(x)_4\}$ next and output the erroneous partition $\{\{bold(x)_1, bold(x)_2, bold(x)_4\}, \{bold(x)_3\}\}$.

  For (2), building $n$ per-cluster heaps $Q_U$ takes $cal(O)(n^2)$ total time via Floyd's linear-time heapify. At merge $t$, updating each of the $n - t - 1$ remaining clusters $C$ performs two deletions and one insertion in $Q_C$ plus one key update in $G$, costing $cal(O)(log n)$ per cluster and $cal(O)((n - t) log n)$ per merge; building $Q_W$ takes $cal(O)(n - t)$. Summing over $t = 1, dots, n - 1$ yields $cal(O)(n^2 log n)$ deterministic worst-case time and $cal(O)(n^2)$ memory.
]

#subsec("E", "Symmetry Invariance, Size-Domination Masking, and Breakdown")

#proposition(title: "Similarity Equivariance, Size-Domination Masking, and 1/n Breakdown", origin: "Original Proof")[
  1. Under any Euclidean similarity $bold(x) |-> a bold(Q) bold(x) + bold(b)$ ($a > 0$, $bold(Q)^top bold(Q) = bold(I)_p$, $bold(b) in bb(R)^p$), UPGMC produces an identical hierarchy with merge heights $tilde(h)_t = a h_t$, while failing invariance under non-isotropic linear maps.
  2. *Size-domination masking*: Merging a cluster $A$ of size $m$ with a cluster $B$ of size $r$ moves the centroid of $A$ by exactly $norm(bold(mu)_(A union B) - bold(mu)_A)_2 = (r / (m + r)) d_C (A, B) = cal(O)(r / m)$ as $m -> infinity$.
  3. *Finite-sample $1/n$ breakdown*: Moving a single point $bold(x)^* in X$ to satisfy $min_(bold(x) in X without \{bold(x)^*\}) norm(bold(x)^* - bold(x))_2 > 2 "diam"(X without \{bold(x)^*\})$ forces $\{bold(x)^*\}$ to remain an isolated singleton until step $t = n - 1$, corrupting every fixed-$k$ cut ($k >= 2$).
]

#proof[
  (1) Linearity of the mean gives $bold(mu)_(a bold(Q) A + bold(b)) = a bold(Q) bold(mu)_A + bold(b)$, and orthogonality gives $norm(a bold(Q) (bold(mu)_A - bold(mu)_B))_2 = a norm(bold(mu)_A - bold(mu)_B)_2$, preserving all distance orderings. Non-isotropic scaling alters Euclidean distances unevenly and can flip merge order.
  (2) Directly from $bold(mu)_(A union B) - bold(mu)_A = (r / (m + r)) (bold(mu)_B - bold(mu)_A)$.
  (3) For any subset $C subset.eq X' := X without \{bold(x)^*\}$, convexity of the mean implies $bold(mu)_C in "conv"(X')$, so $norm(bold(mu)_U - bold(mu)_V)_2 <= "diam"(X')$ for all $U, V subset.eq X'$. By the triangle inequality, for any $bold(x) in X'$, $norm(bold(x)^* - bold(mu)_C)_2 >= norm(bold(x)^* - bold(x))_2 - norm(bold(x) - bold(mu)_C)_2 > 2 "diam"(X') - "diam"(X') = "diam"(X')$. Hence $\{bold(x)^*\}$ can never win a merge against any pair inside $X'$ and remains isolated until $t = n - 1$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

UPGMC occupies a distinctive operational regime in hierarchical clustering: unlike pure dissimilarity methods (Single, Complete, UPGMA, WPGMA), its update law is tied intrinsically to Euclidean geometry (Lemma 1), and its negative coefficient $beta = -n_A n_B / (n_A + n_B)^2$ forbids NN-chain acceleration @mullner2011modern. In practical libraries (such as SciPy `fastcluster`), Müllner's *generic linkage* algorithm maintains a priority queue of candidate nearest-neighbor bounds and triggers queue recomputations only upon encountering inversions; this achieves near-$cal(O)(n^2)$ empirical time when inversions are rare, degrading to $cal(O)(n^3)$ worst-case on inversion-heavy high-dimensional simplices (Theorem 1(3)), whereas the indexed-heap schedule (Algorithm box) guarantees $cal(O)(n^2 log n)$ worst-case time unconditionally.

Memory consumption is $cal(O)(n^2)$ to store the active distance matrix and heap indices (or $cal(O)(n p)$ if centroids are maintained explicitly and pairwise distances are recomputed on the fly in $cal(O)(n^3 p)$ time). Cache locality is poorer than NN-chain WPGMA because heap sift-up/sift-down operations generate random memory accesses across $n$ heaps. Geometrically, UPGMC induces piecewise-spherical, centroid-centered basins: because it removes Ward's size penalty (Theorem 2(2)), large clusters exhibit strong inertia (Proposition 2(2)), steadily absorbing nearby small clusters while leaving isolated outliers unmerged until the top of the tree.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Space], [Euclidean space $bb(R)^p$ (combinatorial on squared $L_2$ matrix $D^((0))$)],
    [Lance--Williams], [$alpha_A = n_A / n_(A B)$, $alpha_B = n_B / n_(A B)$, $beta = -n_A n_B / n_(A B)^2$, $gamma = 0$],
    [Time Complexity], [$cal(O)(n^2 log n)$ worst-case (indexed heaps); $cal(O)(n^3)$ naive],
    [Space Complexity], [$cal(O)(n^2)$ distance matrix + heaps (or $cal(O)(n p)$ explicit centroids)],
    [Reducibility], [Violated ($beta < 0$); inversions possible for all $p >= 2$, $h_(t+1) >= (sqrt(3)/2) h_t$],
    [Ultrametricity], [Not ultrametric in general ($p >= 2$); strictly ultrametric in $p = 1$],
    [Invariance], [Rigid motions and isotropic scaling; non-invariant to affine shear],
    [Inductive Bias], [Spherical barycenter proximity; large-cluster centroid inertia],
  ),
  caption: [Operational and mathematical profile of Centroid Linkage / UPGMC (Alg. 61).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Geospatial Facility Aggregation and Regional Gravity Modeling*: In urban geography and logistics planning, populations or depots are aggregated into regional service centers whose physical center of mass determines transport distance. UPGMC directly merges regions whose mass-weighted geographic barycenters are closest, ensuring that cluster representatives correspond to true physical centers of mass rather than abstract average pairwise distances.
2. *Astronomical Star-Cluster and Halo Merging Trees*: In $N$-body cosmological simulations and astrometry, protostellar clumps or dark-matter subhalos interact primarily through the separation of their centers of mass. When gravitational mass scales with particle count, UPGMC's cardinality-weighted centroid recurrence mirrors exact center-of-mass conservation during hierarchical mergers, naturally capturing the absorption of satellite clumps by massive host halos (Proposition 2(2)).
3. *Vector Quantization and Codebook Tree Construction*: In hierarchical speech and image compression, quantization cells are represented by their Euclidean centroids (reconstructing each cell at minimum mean squared error). UPGMC constructs a bottom-up codebook hierarchy whose merge metric directly tracks inter-codeword Euclidean separation without forcing equal cell populations (unlike Ward's method @ward1963hierarchical), accommodating skewed signal symbol frequencies.
4. *Chemometrics and Spectroscopic Prototype Libraries*: Infrared and NMR spectra of chemical mixtures combine linearly in Euclidean intensity space; averaging replicate spectra yields a physically meaningful denoised mean spectrum (with variance reduced by $1/|A|$). UPGMC groups replicate spectral measurements around true mean prototypes while separating chemically distinct functional groups.
5. *Large-Scale Micro-Cluster Post-Aggregation*: When paired with a streaming front-end such as BIRCH @zhang1996birch that summarizes millions of points into $m << n$ Clustering Feature (CF) micro-clusters storing $(N, bold(L S), S S)$, the unweighted centroid $bold(mu) = bold(L S) / N$ is immediately available in $cal(O)(1)$ time. Running UPGMC over the $m$ micro-clusters reconciles true population weights across heterogeneous micro-cluster sizes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

UPGMC exhibits four critical mathematical failure regimes:
- *Dendrogram inversions and ambiguous tree cuts*: When $h_(t+1) < h_t$ (Theorem 1), cutting the dendrogram at a fixed threshold $tau$ ceases to define a unique nested partition, complicating visual interpretation and tree-based validity indices.
- *Non-convex and hollow manifolds*: For concentric rings, moons, or filaments, a cluster's Euclidean centroid $bold(mu)_A$ can fall in empty space inside another cluster's interior, triggering catastrophic cross-manifold merges.
- *Severe size-domination masking*: By Proposition 2(2), when a large cluster absorbs a smaller nearby cluster, its centroid barely shifts, effectively masking the distinct identity of the smaller subpopulation; Median Linkage (Alg. 62) was explicitly designed to mitigate this effect.
- *Single-outlier partition collapse and non-Euclidean failure*: By Proposition 2(3), a single extreme outlier hijacks an entire cluster slot in any top-$k$ cut, and applying UPGMC to arbitrary non-Euclidean dissimilarities (e.g., Jaccard or edit distance) breaks the Hilbert-space Stewart identity (Lemma 1).

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Physical Center-of-Mass Aggregation], [High: exact mass-weighted barycenter kinematics], [Watch inversions; use Ward @ward1963hierarchical for monotone cuts],
    [Spectral / Prototype Denoising], [High: linear mean prototype reduces noise by $1/|A|$], [Outliers corrupt mean: use PAM / medoids @kaufman1990finding],
    [BIRCH CF Micro-Cluster Merging], [High: exact size-weighted centroid fusion from $(N, bold(L S))$], [Balanced partitions: prefer Ward or $k$-means++ @arthur2007kmeanspp],
    [Non-Convex Manifolds & Moons], [Low: centroids fall off-manifold into empty space], [Use DBSCAN @ester1996dbscan, HDBSCAN @campello2013hdbscan, or Spectral @shi2000ncut],
    [Non-Euclidean / Phylogenetic Matrices], [Low: violates Hilbert-space geometry of Lemma 1], [Use WPGMA (Alg. 60), UPGMA, or Neighbor-Joining @saitou1987neighbor],
  ),
  caption: [Domain suitability and failure-mode matrix for Centroid Linkage (UPGMC).],
)

#sec("V", "Structural Relations and Empirical Behavior")

UPGMC forms the unweighted geometric anchor of the centroid-based Lance--Williams quartet (Ward, UPGMA, UPGMC, WPGMC). The following proposition quantifies its exact asymptotic convergence toward reducibility on imbalanced merges — a sharp structural contrast with Median Linkage (WPGMC, Alg. 62) — and characterizes its relation to $k$-means @lloyd1982least.

#proposition(title: "Asymptotic Reducibility on Imbalanced Merges and Fixed-Point Separation from K-Means", origin: "Original Proof")[
  1. *Imbalanced-merge stability vs. WPGMC*: Suppose UPGMC merges $(A, B)$ with size ratio $r = max(|A|, |B|) / min(|A|, |B|) >= 1$. Then the maximum possible fractional squared-height inversion decays as
     $ (H_t - H_(t+1)) / H_t <= r / (1 + r)^2 <= 1 / r arrow.b 0 quad (r -> infinity), $
     whereas WPGMC (Alg. 62) retains constant $beta = -1/4$ and permits up to $25\%$ squared-height drops at every $r >= 1$.
  2. *Suboptimality relative to $k$-means*: Cutting a UPGMC hierarchy at $k$ clusters does not in general yield a Lloyd stationary partition; however, its merge criterion minimizes $Delta "SSE"(A, B) / (|A| |B| / (|A| + |B|))$, acting as a size-debiased greedy proxy for spherical variance clustering.
]

#proof[
  (1) From Theorem 1(1), $H_(t+1) >= (1 - alpha (1 - alpha)) H_t$ with $alpha = |A| / (|A| + |B|)$. Writing $alpha (1 - alpha) = r / (1 + r)^2$ yields $(H_t - H_(t+1)) / H_t <= r / (1 + r)^2 < 1 / r$. By contrast, WPGMC sets $beta = -1/4$ identically regardless of $r$.
  (2) Let $C_0 subset [-10^(-3), 10^(-3)]$ consist of $10$ distinct points with mean $0$, and set $X = C_0 union \{3.0, 5.0, 8.1\} subset bb(R)$ with $k = 2$. Because all intra-$C_0$ distances are at most $2 dot 10^(-3)$, Steps 1--9 merge $C_0$ to centroid $0$. Step 10 merges $\{3.0, 5.0\}$ ($h_(10) = 2.0 < 3.0 < 3.1$) to centroid $4.0$. Step 11 compares $d_C (C_0, \{3.0, 5.0\}) = 4.0$ against $d_C (\{3.0, 5.0\}, \{8.1\}) = 4.1$, merging $A = C_0 union \{3.0, 5.0\}$ to centroid $mu_A = 8 / 12 = 2 / 3$ and leaving $B = \{8.1\}$ ($mu_B = 8.1$). However, $|5.0 - mu_A| = 13 / 3 approx 4.333 > 3.100 = |5.0 - mu_B|$, so point $5.0$ is strictly closer to the opposing centroid $mu_B$, proving non-stationarity under Lloyd's $k$-means assignment step. The scaling form follows directly from Theorem 2(2).
]

*Empirical Behavior on the ClusterOverview Benchmark*: On the 4-dataset ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), Centroid Linkage (UPGMC) achieves *Iris* $"Sil" = 0.489$, $"ARI" = 0.560$; *Customers* $"Sil" = 0.353$; *Multi-Shape* $"Sil" = 0.452$ ($"ARI" =$ N/A); and *High-D* $"Sil" = 0.234$ ($"ARI" =$ N/A), with a mean runtime of $2.104$ ms. Three quantitative findings mirror our proofs:
- *High Iris Silhouette vs. Lower Iris ARI*: UPGMC attains the highest *Iris* Silhouette (`0.489`) among Algorithms 60--63, yet its *Iris* ARI (`0.560`) trails WPGMC (`0.653`, Alg. 62) by nearly 0.10. This directly illustrates size-domination masking (Proposition 2(2)): once a dominant cluster forms in the overlapping *Versicolor* / *Virginica* region, centroid inertia produces geometrically compact but label-skewed unbalanced clusters that score well on Silhouette while sacrificing ground-truth species alignment.
- *Off-Manifold Degradation on Multi-Shape*: UPGMC drops to $"Sil" = 0.452$ on *Multi-Shape* (identical to WPGMC and well below WPGMA's `0.609`), confirming that centroid compression fails on non-convex shapes whose barycenters lie outside the data manifold.
- *High-D Concentration & Runtime*: On *High-D*, $"Sil" = 0.234$ reflects Euclidean distance concentration @beyer1999nearest and simplex-like inversion pressure (Theorem 1(3)), while the $2.104$ ms runtime aligns closely with WPGMA (`2.209` ms) and WPGMC (`2.225` ms) on moderate-$n$ benchmarks where priority-queue overhead stays modest.
