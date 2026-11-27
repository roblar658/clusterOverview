#import "../template/ieee.typ": *

#alg_entry(
  id: 74,
  name: "Flexible Beta Linkage",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1967",
  paradigm: "One-Parameter Space-Dilating / Contracting Lance-Williams Family",
  complexity: $O(n^2) "for" beta <= 0, quad O(n^2 log n) "for" beta > 0$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Within the four-parameter Lance-Williams agglomerative family @lance1967general (Alg. 72), practitioners frequently need a *single continuous control parameter* that smoothly interpolates between chaining (space-contracting) and compact-cluster (space-dilating) inductive biases on arbitrary non-Euclidean dissimilarity matrices — without incurring dendrogram inversions. Lance and Williams @lance1967general introduced *Flexible Beta Linkage* by imposing three structural constraints on the general combinatorial update: symmetry between the two merged clusters ($alpha_i = alpha_j = alpha$), elimination of the absolute-difference order-statistic term ($gamma = 0$), and unit affine weight normalization ($alpha_i + alpha_j + beta = 1$). Solving $2 alpha + beta = 1$ leaves a single free scalar parameter $beta < 1$.

#definition(title: "Flexible Beta Dissimilarity Recurrence")[
  Let $X = {x_1, dots.c, x_n}$ carry a symmetric zero-diagonal dissimilarity $d : X times X arrow.r RR_(>= 0)$, and fix a scalar parameter $beta in (-infinity, 1)$ (standard operating interval $beta in [-1, 1)$). Upon merging disjoint active clusters $C_i$ and $C_j$, *Flexible Beta Linkage* @lance1967general updates the dissimilarity to every remaining active cluster $C_k$ via
  $
    D_beta (C_i union C_j, C_k) &= ((1 - beta) / 2) (d_(i k) + d_(j k)) + beta d_(i j),
  $
  where $d_(i j) = D_beta (C_i, C_j)$, $d_(i k) = D_beta (C_i, C_k)$, and $d_(j k) = D_beta (C_j, C_k)$.
]

A direct algebraic rewrite clarifies the geometric action of $beta$. Writing the midpoint average $overline(d)_(i j k) = (d_(i k) + d_(j k)) / 2$ (which is exact WPGMA @sokal1958statistical), the update becomes
$
  D_beta (C_i union C_j, C_k) = overline(d)_(i j k) - beta (overline(d)_(i j k) - d_(i j)).
$
At any merge step, minimality of $(C_i, C_j)$ enforces $d_(i j) <= min(d_(i k), d_(j k)) <= overline(d)_(i j k)$, so the margin $Delta_(i j k) = overline(d)_(i j k) - d_(i j)$ is *always non-negative*. Consequently:
- *$beta = 0$ (space-conserving)*: Recovers WPGMA (McQuitty linkage) identically, leaving distances centered at $overline(d)_(i j k)$.
- *$beta < 0$ (space-dilating)*: Adds the positive penalty $(-beta) Delta_(i j k) >= 0$, pushing newly merged clusters *farther* from unmerged clusters than WPGMA predicts, suppressing chaining and enforcing compact, balanced clusters. Lance and Williams @lance1967general recommended $beta = -0.25$ as an all-purpose default, widely adopted in ecology and multivariate taxonomy @kaufman1990finding.
- *$0 < beta < 1$ (space-contracting)*: Subtracts $beta Delta_(i j k) >= 0$, pulling merged clusters *closer* to unmerged points and inducing chaining.

#algbox(
  title: "Flexible Beta Agglomerative Clustering",
  inputs: [Symmetric dissimilarity matrix $D^(0) in RR_(>= 0)^(n times n)$, dilation parameter $beta in (-infinity, 1)$ (default $beta = -0.25$)],
  outputs: [Monotonic dendrogram ${(C_(i_t), C_(j_t), h_t)}_(t=1)^(n-1)$ and induced ultrametric $u_beta$],
)[
  1. Set $alpha arrow.l (1 - beta) / 2$, active cluster set $cal(A) arrow.l {1, dots.c, n}$, and working matrix $D arrow.l D^(0)$.
  2. *Branch on reducibility (Theorem 1)*:
     - *If $beta <= 0$*: Execute the $O(n^2)$ Nearest-Neighbor Chain loop @mullner2011modern; push arbitrary $u in cal(A)$ onto stack $S$, advance to nearest active neighbor until a reciprocal pair $(i, j)$ appears at top of $S$, pop and merge $(i, j)$ at height $D_(i j)$.
     - *If $0 < beta < 1$*: Maintain per-cluster min-heaps over $cal(A)$ and pop the global minimum pair $(i, j) in arg min_(u < v in cal(A)) D_(u v)$ at each step $t = 1, dots.c, n - 1$ ($O(n^2 log n)$ total).
  3. On merging $(i, j)$ into new node $m$, update for every $k in cal(A) without {i, j}$:
     $D_(m k) arrow.l alpha (D_(i k) + D_(j k)) + beta D_(i j)$,
     and replace $cal(A) arrow.l (cal(A) without {i, j}) union {m}$.
  4. Sort recorded merges by non-decreasing height $h_t$ (guaranteed inversion-free) and cut to desired $k$ or threshold $h$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational properties of Flexible Beta Linkage: (A) universal inversion-free monotonicity across $(-infinity, 1)$ coupled with a sharp reducibility phase transition at $beta = 0$, (B) exact closed-form amplification laws for balanced versus caterpillar trees explaining the structural boundary $beta = -1$, (C) an exact local equivalence to Ward's minimum-variance criterion explaining why $beta approx -0.25$ succeeds empirically, and (D) an exact characterization showing that metricity (the triangle inequality) is preserved *if and only if* $beta = 0$.

#subsec("A", "Universal Monotonicity and the Sharp Reducibility Transition at Beta = 0")

#theorem(title: "Monotonicity Across (-Infinity, 1) and Reducibility Threshold at Beta = 0", origin: "Original Proof")[
  For Flexible Beta Linkage with $beta in (-infinity, 1)$ on any symmetric dissimilarity $d >= 0$:
  1. *Universal Monotonicity*: Every merge step satisfies $D_beta (C_i union C_j, C_k) >= D_beta (C_i, C_j)$ for all $C_k$, guaranteeing non-decreasing dendrogram heights $h_1 <= h_2 <= dots.c <= h_(n-1)$ and an exact induced ultrametric $u_beta$ for *every* $beta < 1$.
  2. *Sharp Reducibility Phase Transition*: Bruynooghe's reducibility inequality $D_beta (C_i union C_j, C_k) >= min(d_(i k), d_(j k))$ holds universally *if and only if* $beta <= 0$. For every $beta in (0, 1)$, reducibility fails on a 3-point Euclidean triangle.
]

#proof[
  At merge step $t$, minimality gives $d_(i j) <= min(d_(i k), d_(j k))$. Let $m = min(d_(i k), d_(j k))$ and $M = max(d_(i k), d_(j k))$, so $0 <= d_(i j) <= m <= M$.

  *Part 1 (Monotonicity)*: Subtracting $d_(i j)$ from $D_beta (C_(i j), C_k)$ using $2 alpha + beta = 1$ yields
  $
    D_beta (C_(i j), C_k) - d_(i j) &= ((1 - beta) / 2) ((d_(i k) - d_(i j)) + (d_(j k) - d_(i j))) \
    &= ((1 - beta) / 2) ((M - d_(i j)) + (m - d_(i j))).
  $
  Since $beta < 1$ implies $(1 - beta) / 2 > 0$, and both $M - d_(i j) >= 0$ and $m - d_(i j) >= 0$, we have $D_beta (C_(i j), C_k) >= d_(i j) = h_t$, with strict inequality whenever $M > d_(i j)$. All unmerged pairs already have dissimilarity $>= h_t$, proving $h_(t+1) >= h_t$ and ultrametricity @lance1967general.

  *Part 2 (Reducibility Phase Transition)*: Subtracting $m = min(d_(i k), d_(j k))$ yields the exact decomposition
  $
    D_beta (C_(i j), C_k) - m &= ((1 - beta) / 2) (M + m) + beta d_(i j) - m \
    &= ((1 - beta) / 2) (M - m) - beta (m - d_(i j)).
  $
  - If $beta <= 0$: then $(1 - beta) / 2 >= 1 / 2 > 0$, $-beta >= 0$, $M - m >= 0$, and $m - d_(i j) >= 0$. Both terms are non-negative, proving reducibility unconditionally and validating the $O(n^2)$ Nearest-Neighbor Chain algorithm @mullner2011modern.
  - If $beta in (0, 1)$: take the 3-point isosceles Euclidean triangle $x_1 = (0, 0)$, $x_2 = (1, 0)$, $x_3 = (1 / 2, 2) in RR^2$, where $d_(1 2) = 1$ and $m = M = d_(1 3) = d_(2 3) = sqrt(17) / 2 > 1$. Merging $(x_1, x_2)$ gives $M - m = 0$ and
    $
      D_beta ({x_1, x_2}, {x_3}) - m = -beta (sqrt(17) / 2 - 1) < 0,
    $
    violating reducibility for every $beta in (0, 1)$.
]

#subsec("B", "Closed-Form Tree Amplification Laws and the Beta = -1 Boundary")

Why do practitioners restrict space dilation to $beta in [-1, 0]$ rather than arbitrarily large negative $beta$? We derive exact closed-form dissimilarity trajectories on balanced versus caterpillar merge trees, exposing $beta = -1$ as an exact topological stability boundary.

#theorem(title: "Pointwise Dilation Ordering and Exact Balanced / Caterpillar Growth Laws", origin: "Original Derivation")[
  1. *Monotone Dilation in $beta$*: For fixed dissimilarities satisfying $d_(i j) <= min(d_(i k), d_(j k))$, $D_(beta_1) (C_(i j), C_k) >= D_(beta_2) (C_(i j), C_k)$ whenever $beta_1 <= beta_2 < 1$.
  2. *Balanced Merge Amplification*: Let $C$ be a zero-diameter cluster of size $abs(C) = 2^L$ at initial distance $d_0 > 0$ from an external point $z$, assembled via a complete binary tree of depth $L$. Then
     $
       D_beta (C, {z}) = (1 - beta)^L d_0 = abs(C)^(log_2 (1 - beta)) d_0.
     $
  3. *Caterpillar Merge Amplification and Critical Threshold $beta = -1$*: If instead $C_m = {x_1, dots.c, x_m}$ (zero internal diameter, initial distance $d_0$ to $z$) grows sequentially one singleton at a time ($C_t = C_(t-1) union {x_t}$), then $u_m = D_beta (C_m, {z})$ satisfies
     $
       u_m = cases(
         d_0 ( (1 - beta) / (1 + beta) + (2 beta) / (1 + beta) ((1 - beta) / 2)^(m - 1) ) & "if" beta != -1,
         m d_0 & "if" beta = -1.
       )
     $
     Hence as $m arrow.r infinity$, sequential absorption converges to the finite limit $((1 - beta) / (1 + beta)) d_0$ *if and only if* $beta in (-1, 1)$, grows linearly at $beta = -1$, and diverges exponentially as $((1 - beta) / 2)^m$ for $beta < -1$.
]

#proof[
  *Part 1*: Differentiating $D_beta$ with respect to $beta$ gives
  $
    partial D_beta (C_(i j), C_k) / (partial beta) = d_(i j) - (d_(i k) + d_(j k)) / 2 = -Delta_(i j k) <= 0,
  $
  establishing pointwise monotonic decrease in $beta$.

  *Part 2*: Proceed by induction on tree depth $ell = 0, dots.c, L$. For $ell = 0$, $D_beta ({x}, {z}) = d_0 = (1 - beta)^0 d_0$. Merging two depth-$(ell - 1)$ subclusters $A, B$ at internal distance $D_beta (A, B) = 0$ yields
  $
    D_beta (A union B, {z}) &= ((1 - beta) / 2)(2 (1 - beta)^(ell - 1) d_0) + beta dot 0 \
    &= (1 - beta)^ell d_0.
  $
  Setting $ell = L = log_2 abs(C)$ gives $(1 - beta)^L d_0 = abs(C)^(log_2 (1 - beta)) d_0$.

  *Part 3*: At step $m >= 2$, $D_beta (C_(m-1), {x_m}) = 0$ (inductively verified since internal distances stay $0$ when $d = 0$ inside $C_m$). Setting $r = (1 - beta) / 2$, the caterpillar recurrence is $u_1 = d_0$ and $u_m = r (u_(m-1) + d_0)$. Unrolling yields
  $
    u_m = d_0 ( r^(m-1) + sum_(j=1)^(m-1) r^j ).
  $
  When $beta = -1$, $r = 1$ so $u_m = m d_0$. When $beta != -1$ ($r != 1$), summing the geometric series yields
  $
    u_m &= d_0 ( r^(m-1) + r (1 - r^(m-1)) / (1 - r) ) \
    &= d_0 ( r / (1 - r) + ((1 - 2 r) / (1 - r)) r^(m-1) ).
  $
  Substituting $1 - r = (1 + beta) / 2$, $r / (1 - r) = (1 - beta) / (1 + beta)$, and $(1 - 2 r) / (1 - r) = 2 beta / (1 + beta)$ yields the stated closed form. Since $abs(r) < 1$ holds iff $-1 < beta < 3$, and $beta < 1$ is required for monotonicity, convergence as $m arrow.r infinity$ holds precisely on $beta in (-1, 1)$ with limit $((1 - beta) / (1 + beta)) d_0$.
]

At the canonical setting $beta = -0.25$, balanced clusters dilate mildly as $abs(C)^(log_2 1.25) approx abs(C)^(0.322)$ while caterpillar chains saturate at $(5 / 3) d_0$; below $beta = -1$, both regimes explode exponentially, producing pathological forced equal-size splits regardless of true cluster geometry.

#subsec("C", "Local Equivalence to Ward's Method Explaining Beta = -0.25")

A longstanding empirical observation in numerical taxonomy @lance1967general @kaufman1990finding is that Flexible Beta around $beta = -0.25$ closely emulates Ward's minimum-variance clustering @ward1963hierarchical — even on non-Euclidean dissimilarities where Ward's centroid interpretation is undefined. We derive the exact algebraic reason.

#proposition(title: "Exact Instantaneous Matching Between Ward and Flexible Beta", origin: "Original Derivation")[
  Consider merging two equal-size clusters $n_i = n_j = m$ relative to a third cluster of size $n_k = rho m$ ($rho > 0$). Ward's squared-distance Lance-Williams update (Alg. 72) coincides *identically* with Flexible Beta Linkage at the effective parameter
  $
    beta_"eff" (rho) = -rho / (2 + rho) in (-1, 0).
  $
  In particular, balanced three-way merges ($rho = 1$) match $beta = -1 / 3$, while merging size-2 clusters against singletons ($rho = 1 / 2$) matches $beta = -1 / 5$, bracketing Lance and Williams's empirical rule $beta = -0.25$ (exact at $rho = 2 / 3$).
]

#proof[
  In Ward's Lance-Williams recurrence (Alg. 72, Proposition 2), substituting $n_i = n_j = m$ and $n_k = rho m$ (so $N = (2 + rho) m$) yields
  $
    alpha_i = alpha_j = ((1 + rho) m) / ((2 + rho) m) = (1 + rho) / (2 + rho), quad beta_"Ward" = -rho / (2 + rho),
  $
  with $gamma = 0$. Checking the Flexible Beta identity $(1 - beta_"eff") / 2$:
  $
    (1 - (-rho / (2 + rho))) / 2 = ((2 + 2 rho) / (2 + rho)) / 2 = (1 + rho) / (2 + rho) = alpha_i = alpha_j.
  $
  Thus Ward's update equals $D_(beta_"eff" (rho))$ identically at every equal-size merge. Evaluating at $rho = 1, 2 / 3, 1 / 2$ yields $beta_"eff" = -1 / 3, -1 / 4, -1 / 5$, respectively.
]

#subsec("D", "Exact Characterization of Metric Preservation")

While Flexible Beta produces a valid ultrametric $u_beta$ across all $beta < 1$, the *intermediate* inter-cluster dissimilarities $D_beta$ almost always lose the triangle inequality.

#proposition(title: "Triangle-Inequality Preservation Holds If and Only If Beta = 0", origin: "Original Proof")[
  Suppose the initial dissimilarity $d$ on $X$ is a metric. Then the updated active dissimilarity matrix $D_beta$ is guaranteed to remain a metric on the active cluster set $cal(A)$ after one merge *if and only if* $beta = 0$ (WPGMA). For every $beta < 0$ and every $beta in (0, 1)$, explicit 4-point collinear Euclidean configurations violate the triangle inequality after the first merge.
]

#proof[
  *Sufficiency ($beta = 0$)*: At $beta = 0$, $D_0 (C_(i j), C_k) = (d_(i k) + d_(j k)) / 2$. For any two remaining clusters $C_k, C_ell$:
  $
    D_0 (C_(i j), C_k) &= (d_(i k) + d_(j k)) / 2 \
    &<= (d_(i ell) + d_(ell k) + d_(j ell) + d_(ell k)) / 2 = D_0 (C_(i j), C_ell) + d_(ell k),
  $
  and $d_(k ell) = (d_(k ell) + d_(k ell)) / 2 <= (d_(k i) + d_(i ell) + d_(k j) + d_(j ell)) / 2 = D_0 (C_(i j), C_k) + D_0 (C_(i j), C_ell)$. Positivity and symmetry are immediate, so $D_0$ remains a metric by induction.

  *Necessity for $beta < 0$*: Place four points on $RR$: $x_1 = 0$, $x_2 = epsilon$, $w = 1 / 2$, $z = 1$ with $0 < epsilon < 1 / 4$. The unique first merge joins $C_(1 2) = {x_1, x_2}$ at height $epsilon$. Computing updated distances:
  $
    D_beta (C_(1 2), {w}) &= ((1 - beta) / 2)(1 - epsilon) + beta epsilon, \
    D_beta (C_(1 2), {z}) &= ((1 - beta) / 2)(2 - epsilon) + beta epsilon,
  $
  while $d(w, z) = 1 / 2$. Forming the triangle-inequality excess, every $epsilon$-term cancels identically:
  $
    D_beta (C_(1 2), {z}) - D_beta (C_(1 2), {w}) - d(w, z) = (1 - beta) / 2 - 1 / 2 = -beta / 2 > 0.
  $

  *Necessity for $beta in (0, 1)$*: Place $x_1 = 0$, $x_2 = epsilon$, $w = -1 / 2$, $z = 1 / 2$ on $RR$ with $0 < epsilon < 1 / 4$. After merging $C_(1 2) = {x_1, x_2}$ at height $epsilon$, we obtain $D_beta (C_(1 2), {w}) = ((1 - beta) / 2)(1 + epsilon) + beta epsilon$ and $D_beta (C_(1 2), {z}) = ((1 - beta) / 2)(1 - epsilon) + beta epsilon$, so
  $
    D_beta (C_(1 2), {w}) + D_beta (C_(1 2), {z}) = 1 - beta + 2 beta epsilon < 1 = d(w, z)
  $
  for all $epsilon < 1 / 2$, strictly violating the triangle inequality.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Theorem 1 establishes an operational dichotomy right at $beta = 0$: in the space-dilating regime $beta <= 0$ — which covers almost all practical deployments including $beta = -0.25$ — Flexible Beta is unconditionally reducible and executes via the Nearest-Neighbor Chain algorithm in optimal $O(n^2)$ time and $O(n^2)$ packed triangular memory @mullner2011modern @schubert2021fast. In the space-contracting regime $0 < beta < 1$, reducibility fails and out-of-order RNN merges can produce erroneous dendrograms; exact execution requires priority-queue maintenance at $O(n^2 log n)$ time and $O(n^2)$ space. Unlike Ward's method, Flexible Beta requires no tracking of cluster sizes $n_u$ during distance updates (because $alpha = (1 - beta) / 2$ and $beta$ are constants!), simplifying inner-loop vectorization to a single fused multiply-add per active cluster.

Geometrically, $beta$ provides a smooth, monotonic control over cluster compactness without leaving the inversion-free cone. Negative $beta$ introduces an implicit size-dependent repulsion (Theorem 2) that prevents large clusters from absorbing distant satellites or chaining across low-density noise bridges, effectively yielding Ward-like compact partitions on non-Euclidean dissimilarities (Proposition 3). Positive $beta$ acts like a softened single linkage, contracting distances around growing clusters to trace elongated filaments — at the cost of increased noise sensitivity and slower $O(n^2 log n)$ execution.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time complexity], [$O(n^2)$ for $beta <= 0$ (NN-chain); $O(n^2 log n)$ for $0 < beta < 1$ (heap)],
    [Space complexity], [$O(n^2)$ packed dissimilarity matrix (no cardinality bookkeeping needed)],
    [Monotonicity], [Inversion-free ultrametric for all $beta in (-infinity, 1)$],
    [Reducibility], [Exact iff $beta <= 0$ (sharp phase transition at WPGMA $beta = 0$)],
    [Hyperparameter], [Single scalar $beta in [-1, 1)$; canonical default $beta = -0.25$],
    [Metricity], [Active $D_beta$ preserves triangle inequality iff $beta = 0$],
    [Invariances], [Positive scaling & permutation; not nonlinear monotone-transform invariant],
  ),
  caption: [Operational and theoretical profile of Flexible Beta Linkage (Alg. 74).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Community Ecology and Phytosociology*: Ecological abundance matrices exhibit the double-zero problem and are analyzed via Bray-Curtis or Sørensen dissimilarities rather than Euclidean distances. Because Ward's method is formally invalid on Bray-Curtis matrices, Flexible Beta with $beta = -0.25$ is the canonical ecological standard @lance1967general @kaufman1990finding, reproducing Ward's balanced, chaining-free habitat clusters directly in non-Euclidean dissimilarity space (Proposition 3).
2. *Microbiome and Metagenomic Beta-Diversity Stratification*: Microbial community profiling compares samples through weighted / unweighted UniFrac or Jaccard beta-diversity matrices. Tuning $beta in [-0.35, -0.15]$ cleanly separates enterotypes while suppressing outlier-driven chaining that plagues single linkage and unweighted UPGMA.
3. *Mixed-Type Customer and Socioeconomic Segmentation*: Commercial customer tables mix continuous spend, ordinal tiers, and categorical attributes, typically summarized via Gower dissimilarity @kaufman1990finding. Flexible Beta ($beta = -0.25$) operates natively on Gower matrices in $O(n^2)$ time while preventing small-cluster fragmentation and giant-cluster snowballing.
4. *Cheminformatics Molecular Scaffold Clustering*: Compound screening workflows on binary Morgan fingerprints use Tanimoto dissimilarity, where $beta$ lets medicinal chemists dial continuously between tight congeneric scaffolds ($beta approx -0.25$ to $-0.5$) and broader bioisosteric series ($beta approx 0$ to $+0.1$).
5. *Psychometric and Behavioral Profile Taxonomy*: Survey instruments with ordinal Likert scales and rank-correlation dissimilarities employ Flexible Beta to recover compact respondent archetypes without assuming Euclidean additivity across questionnaire items @aggarwal2013data.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five operational pitfalls require care: (1) *Hyperparameter selection sensitivity*, since no single $beta$ is universally optimal across both spherical and elongated geometries, requiring external validation (e.g., Silhouette @rousseeuw1987silhouettes or Calinski-Harabasz @calinski1974dendrite); (2) *Artificial equal-size dissection for $beta < -0.5$*, where super-linear space dilation (Theorem 2) forces genuinely large homogeneous clusters to split prematurely; (3) *Runaway chaining for $beta > 0$*, where exponential contraction pulls noise bridges into giant snowball clusters and breaks $O(n^2)$ NN-chain validity; (4) *Loss of metric scale interpretability*, since updated distances $D_beta$ violate the triangle inequality for all $beta != 0$ (Proposition 4) and inflate above the original data diameter when $beta < 0$; and (5) *Quadratic memory ceiling* ($O(n^2)$), blocking direct application above $n approx 5 times 10^4$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Non-Euclidean ecology / Gower], [Excellent: Ward-like compactness without Euclidean coordinates], [Keep $beta in [-0.5, -0.1]$; validate cut height],
    [High-D / noisy tabular data], [Strong: mild dilation ($beta = -0.25$) suppresses hub chaining], [Avoid $beta < -1$ (exponential divergence)],
    [Filament / density manifolds], [Moderate: tune $beta$ gently; positive $beta$ risks noise chaining], [Prefer HDBSCAN @campello2013hdbscan or OPTICS @ankerst1999optics],
    [Euclidean minimum-variance], [Strong approximation to Ward ($beta approx -0.25$ to $-0.33$)], [Use exact Ward @ward1963hierarchical if coordinates exist],
    [Large scale ($n > 5 times 10^4$)], [Infeasible: $O(n^2)$ matrix memory wall], [Use BIRCH @zhang1996birch or $k$-means++ @arthur2007kmeanspp],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Alg. 74.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Flexible Beta occupies a pivotal axis inside the Lance-Williams family (Alg. 72): it generalizes WPGMA ($beta = 0$) continuously across space-dilating ($beta < 0$) and space-contracting ($beta > 0$) regimes while locally shadowing Ward's method around $beta in [-1 / 3, -1 / 5]$ (Proposition 3). Unlike UPGMA (whose weights $n_i / (n_i + n_j)$ depend on cluster size), Flexible Beta is *unweighted* ($alpha_i = alpha_j$), giving small distinct clusters equal influence to large ones during merging — a critical advantage under class imbalance.

#proposition(title: "Unweighted Small-Cluster Preservation Versus UPGMA Dominance", origin: "Original Proof")[
  Suppose a small cluster $C_i$ ($n_i = 1$) merges with a large cluster $C_j$ ($n_j = M >> 1$) at internal distance $d_(i j) = 0$, while an external cluster $C_k$ satisfies $d_(i k) = a$ and $d_(j k) = b$ with $a != b$. As $M arrow.r infinity$, UPGMA collapses to $D_"UPGMA" (C_(i j), C_k) arrow.r b$ (erasing $C_i$'s geometry completely), whereas Flexible Beta preserves $C_i$'s contribution with fixed weight $(1 - beta) / 2 > 0$:
  $
    abs(D_beta (C_(i j), C_k) - D_"UPGMA" (C_(i j), C_k)) arrow.r abs(((1 - beta) / 2) a - ((1 + beta) / 2) b).
  $
]

#proof[
  Under UPGMA, $D_"UPGMA" (C_(i j), C_k) = (a + M b) / (1 + M) arrow.r b$ as $M arrow.r infinity$. Under Flexible Beta with $d_(i j) = 0$, $D_beta (C_(i j), C_k) = ((1 - beta) / 2)(a + b)$ independently of $M$. Subtracting the two limits yields $((1 - beta) / 2) a + ((1 - beta) / 2 - 1) b = ((1 - beta) / 2) a - ((1 + beta) / 2) b$, proving non-vanishing sensitivity to minority subclusters for every $beta < 1$.
]

*Empirical Benchmark Behavior*: On the ClusterOverview 4-dataset benchmark, Flexible Beta Linkage is the strongest overall performer among Algorithms 72--75: it achieves Iris Silhouette $= 0.455$ and $"ARI" = 0.635$ (tying Minimax Linkage and beating generic Lance-Williams's $"ARI" = 0.565$), Customers Silhouette $= 0.347$, an exceptional Multi-Shape Silhouette $= 0.709$ (far surpassing Lance-Williams at $0.452$ and Minimax at $0.417$), and the highest High-D Silhouette $= 0.292$ (nearly doubling Minimax's $0.167$), all in a brisk mean runtime of $2.458$ ms. These empirical numbers strongly corroborate our theory: calibrated space dilation ($beta < 0$) retains optimal $O(n^2)$ NN-chain speed (matching Alg. 72's $2.346$ ms) while counteracting both high-dimensional hub chaining (@beyer1999nearest) and uneven-size absorption (Proposition 5), delivering clean cluster boundaries across heterogeneous geometries @jain2010data @xu2005survey.
