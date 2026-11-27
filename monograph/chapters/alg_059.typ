#import "../template/ieee.typ": *

#alg_entry(
  id: 59,
  name: "Ward's Minimum-Variance Method",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1963",
  paradigm: "Agglomerative Variance-Minimization (ESS Greedy)",
  complexity: $O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Proposed by Joe H. Ward Jr. @ward1963hierarchical, *Ward's Method* formulates hierarchical agglomerative clustering not through an ad-hoc inter-cluster distance heuristic, but via greedy step-by-step minimization of a global statistical loss function: the *Error Sum of Squares* (ESS, equivalently within-cluster sum of squared deviations from cluster centroids). Given a Euclidean dataset $X = {x_1, dots, x_n} subset bb(R)^p$ and a partition $cal(C) = {C_1, dots, C_k}$, define the centroid $mu_C = (1 / abs(C)) sum_(x in C) x$ and the partition objective
$
  "ESS"(cal(C)) = sum_(j=1)^k "ESS"(C_j), quad "ESS"(C) = sum_(x in C) norm(x - mu_C)_2^2.
$
Starting from $cal(C)_n = {{x_1}, dots, {x_n}}$ where $"ESS"(cal(C)_n) = 0$, every merge of two disjoint active clusters $A, B in cal(C)$ strictly increases (or preserves) total within-cluster variance by the exact nonnegative increment
$
  Delta(A, B) &= "ESS"(A union B) - "ESS"(A) - "ESS"(B) \
  &= (abs(A) abs(B)) / (abs(A) + abs(B)) norm(mu_A - mu_B)_2^2.
$
At each step $t = 1, dots, n - 1$, Ward's method greedily selects the pair $(A_t, B_t) in arg min_(A != B) Delta(A, B)$ minimizing the immediate ESS increase. To align the height scale with Euclidean distance on singletons (where $abs(A) = abs(B) = 1$ gives $Delta({x}, {y}) = (1 / 2) norm(x - y)_2^2$), modern implementations @mullner2011modern define the *Ward dissimilarity* $d_W (A, B) = sqrt(2 Delta(A, B))$, so that $d_W ({x}, {y}) = norm(x - y)_2$ and $D_W (A, B) = d_W (A, B)^2 = 2 Delta(A, B)$ satisfies the Lance--Williams recurrence @lance1967general:
$
  D_W (A union B, C) &= (n_A + n_C) / N D_W (A, C) \
  &quad + (n_B + n_C) / N D_W (B, C) - n_C / N D_W (A, B),
$
where $n_A = abs(A)$, $n_B = abs(B)$, $n_C = abs(C)$, and $N = n_A + n_B + n_C$. Notice the *negative* coefficient $beta = -n_C / N < 0$: remarkably, despite this negative term (which causes inversions in Centroid and Median linkage), Ward's method remains strictly reducible and inversion-free.

#definition(title: "Ward's ESS Hierarchy and Ward Ultrametric")[
  Given $X = {x_1, dots, x_n} subset bb(R)^p$, *Ward's hierarchy* is the nested sequence of partitions $cal(C)_n, dots, cal(C)_1$ obtained by merging $(A_t, B_t) in arg min_(A != B in cal(C)_(n-t+1)) Delta(A, B)$ with variance increment $Delta_t = Delta(A_t, B_t)$ and Euclidean-scale merge height $h_t = sqrt(2 Delta_t)$. Its induced cophenetic dissimilarity $u_W (x, y) = h_(t(x, y))$ is the merge height uniting $x$ and $y$.
]

#algbox(
  title: "Exact O(n^2) Ward Clustering via Nearest-Neighbor Chain",
  inputs: [Points $X = {x_1, dots, x_n} subset bb(R)^p$ (or squared Euclidean matrix $D^(2)$); tie-breaking order $prec$; target $k$.],
  outputs: [Exact Ward merge sequence ${(A_t, B_t, h_t)}_(t=1)^(n-1)$ and $k$-partition $cal(C)_k$.],
)[
  1. Initialize $cal(C) <- {{x_1}, dots, {x_n}}$, sizes $n_({x_i}) <- 1$, $D_W ({x_i}, {x_j}) <- norm(x_i - x_j)_2^2$, stack $S <- []$, $t <- 1$.
  2. *While* $abs(cal(C)) > 1$:
     - If $S$ is empty, push any $C in cal(C)$ onto $S$.
     - Let $A$ be the top of $S$; find $B = arg min_(C in cal(C) without {A}) D_W (A, C)$ (tie-breaking via $prec$).
     - *If* $B$ is on $S$ (forming an RNN pair immediately below $A$):
       - Pop $A, B$ from $S$; record merge $(A_t, B_t, h_t) <- (A, B, sqrt(D_W (A, B)))$.
       - Form $U <- A union B$ with $n_U <- n_A + n_B$ and update $cal(C) <- (cal(C) without {A, B}) union {U}$.
       - For every $C in cal(C) without {U}$, set $N <- n_U + n_C$ and update
         $D_W (U, C) <- ((n_A + n_C) D_W (A, C) + (n_B + n_C) D_W (B, C) - n_C D_W (A, B)) / N$.
       - Increment $t <- t + 1$.
     - *Else:* push $B$ onto $S$.
  3. Sort merges by ascending $h_t$ (guaranteed non-decreasing) and execute the first $n - k$ merges to output $cal(C)_k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Exact Variance Decomposition and Lance--Williams Recurrence")

We derive from first principles Ward's merge increment formula, total ANOVA scatter conservation, and the exact centroid-free Lance--Williams recurrence.

#theorem(title: "Parallel-Axis Variance Identity and Exact Lance--Williams Derivation", origin: "Original Derivation")[
  Let $A, B, C subset bb(R)^p$ be finite disjoint point sets with sizes $n_A, n_B, n_C$ and centroids $mu_A, mu_B, mu_C$, and set $N = n_A + n_B + n_C$. Then:
  1. *Merge increment:* $Delta(A, B) = (n_A n_B / (n_A + n_B)) norm(mu_A - mu_B)_2^2$.
  2. *ANOVA telescoping:* For any complete agglomerative sequence with increments $Delta_1, dots, Delta_(n-1)$, the $k$-cluster partition $cal(C)_k$ satisfies $"ESS"(cal(C)_k) = sum_(t=1)^(n-k) Delta_t$, and total scatter obeys $"ESS"({X}) = sum_(t=1)^(n-1) Delta_t$.
  3. *Exact Lance--Williams recurrence:*
     $
       Delta(A union B, C) &= (n_A + n_C) / N Delta(A, C) + (n_B + n_C) / N Delta(B, C) \
       &quad - n_C / N Delta(A, B).
     $
]

#proof[
  *Part 1:* Let $U = A union B$ and $n_U = n_A + n_B$. The merged centroid is $mu_U = (n_A mu_A + n_B mu_B) / n_U$, whence $mu_A - mu_U = (n_B / n_U)(mu_A - mu_B)$ and $mu_B - mu_U = -(n_A / n_U)(mu_A - mu_B)$. For every $x in A$, write $x - mu_U = (x - mu_A) + (mu_A - mu_U)$. Summing squared norms over $A$ and using $sum_(x in A) (x - mu_A) = 0$ eliminates the cross term (parallel-axis theorem):
  $
    sum_(x in A) norm(x - mu_U)_2^2 = "ESS"(A) + n_A norm(mu_A - mu_U)_2^2.
  $
  Adding the symmetric identity over $B$ yields
  $
    Delta(A, B) &= (n_A n_B^2 + n_B n_A^2) / n_U^2 norm(mu_A - mu_B)_2^2 \
    &= (n_A n_B) / (n_A + n_B) norm(mu_A - mu_B)_2^2.
  $

  *Part 2:* Since $"ESS"(cal(C)_(n-t)) - "ESS"(cal(C)_(n-t+1)) = Delta_t$ and $"ESS"(cal(C)_n) = 0$, telescoping over $t = 1, dots, n - k$ immediately yields $"ESS"(cal(C)_k) = sum_(t=1)^(n-k) Delta_t$, and setting $k = 1$ yields $"ESS"({X}) = sum_(t=1)^(n-1) Delta_t$.

  *Part 3:* Write $mu_U - mu_C = (n_A / n_U)(mu_A - mu_C) + (n_B / n_U)(mu_B - mu_C)$. Using the exact Euclidean polarization identity $norm(alpha u + (1 - alpha) v)_2^2 = alpha norm(u)_2^2 + (1 - alpha) norm(v)_2^2 - alpha (1 - alpha) norm(u - v)_2^2$ with $alpha = n_A / n_U$, $u = mu_A - mu_C$, $v = mu_B - mu_C$ (so $u - v = mu_A - mu_B$), we obtain Stewart's centroid identity:
  $
    norm(mu_U - mu_C)_2^2 &= n_A / n_U norm(mu_A - mu_C)_2^2 \
    &quad + n_B / n_U norm(mu_B - mu_C)_2^2 - (n_A n_B) / n_U^2 norm(mu_A - mu_B)_2^2.
  $
  Multiplying both sides by $(n_U n_C) / N$ gives
  $
    Delta(U, C) &= (n_A n_C) / N norm(mu_A - mu_C)_2^2 + (n_B n_C) / N norm(mu_B - mu_C)_2^2 \
    &quad - n_C / N (n_A n_B) / n_U norm(mu_A - mu_B)_2^2.
  $
  Finally, substituting $n_A n_C norm(mu_A - mu_C)_2^2 = (n_A + n_C) Delta(A, C)$, $n_B n_C norm(mu_B - mu_C)_2^2 = (n_B + n_C) Delta(B, C)$, and $(n_A n_B / n_U) norm(mu_A - mu_B)_2^2 = Delta(A, B)$ yields the exact Lance--Williams recurrence. Multiplying by $2$ establishes the recurrence for $D_W$.
]

#subsec("B", "Ward Reducibility, Non-Inversion, and NN-Chain Exactness")

Centroid linkage (UPGMC) also features a negative Lance--Williams coefficient and notoriously suffers from *dendrogram inversions* ($h_(t+1) < h_t$), invalidating ultrametricity and breaking NN-Chain. We prove why Ward's specific negative coefficient preserves reducibility exactly.

#theorem(title: "Ward Reducibility and Strict Non-Inversion", origin: "Original Proof")[
  For any disjoint $A, B, C subset bb(R)^p$, whenever $Delta(A, B) <= min { Delta(A, C), Delta(B, C) }$, we have
  $
    Delta(A union B, C) >= min { Delta(A, C), \, Delta(B, C) }.
  $
  Consequently: (i) Ward merge heights are monotonically non-decreasing ($Delta_1 <= dots <= Delta_(n-1)$ and $h_1 <= dots <= h_(n-1)$); (ii) the cophenetic metric $u_W$ is always a valid ultrametric on $X$; and (iii) the NN-Chain algorithm computes the exact Ward dendrogram in $Theta(n^2)$ time @mullner2011modern.
]

#proof[
  Let $m = min { Delta(A, C), Delta(B, C) }$. By assumption, $Delta(A, C) >= m$, $Delta(B, C) >= m$, and $-Delta(A, B) >= -m$. Substituting these three inequalities into Theorem 1 (Part 3), noting that $(n_A + n_C) / N > 0$, $(n_B + n_C) / N > 0$, and $n_C / N > 0$, gives
  $
    Delta(A union B, C) &>= (n_A + n_C) / N m + (n_B + n_C) / N m - n_C / N m \
    &= (n_A + n_B + n_C) / N m = m.
  $
  The exact cancellation $(n_A + n_C) + (n_B + n_C) - n_C = N$ establishes Bruynooghe's Reducibility unconditionally. (By contrast, Centroid linkage's coefficients sum to $1 - n_A n_B / (n_A + n_B)^2 < 1$, breaking the lower bound.)

  By induction on step $t$, at step $t + 1$ every existing pair $(P, Q)$ satisfies $Delta(P, Q) >= Delta_t$, and every newly created pair $(A_t union B_t, C)$ satisfies $Delta(A_t union B_t, C) >= min { Delta(A_t, C), Delta(B_t, C) } >= Delta_t$. Hence $Delta_(t+1) >= Delta_t$, guaranteeing zero inversions, ultrametricity of $u_W$, and invariance of reciprocal nearest-neighbor pairs under disjoint merges, so NN-Chain returns the exact Ward hierarchy in $Theta(n^2)$ time.
]

#subsec("C", "Size-Balancing Bias and Isotropic Gaussian MLE Equivalence")

Ward's strong empirical performance on well-formed Euclidean clusters stems directly from two structural properties: a concave harmonic cardinality penalty and an exact Gaussian likelihood equivalence.

#proposition(title: "Harmonic Size-Balance Bias and Spherical Gaussian MLE Step", origin: "Original Derivation")[
  1. *Concave size-balancing bias:* For fixed centroid separation $norm(mu_A - mu_B)_2 = R > 0$ and fixed combined cardinality $n_A + n_B = M$, the merge cost $Delta(A, B) = (n_A (M - n_A) / M) R^2$ is strictly concave in $n_A$, minimized at $n_A = 1$ ($Delta approx R^2$) and maximized at equal sizes $n_A = n_B = M / 2$ ($Delta = (M / 4) R^2$). Hence Ward resists merging two large clusters at upper levels, biasing $k$-partitions toward balanced cluster cardinalities.
  2. *Homoscedastic Gaussian MLE equivalence:* Under the fixed-partition spherical Gaussian model $x_i tilde.op cal(N)(mu_(c(i)), sigma^2 I_p)$ @fraley2002mclust, merging $A$ and $B$ changes the profile log-likelihood by
     $
       ell(cal(C)') - ell(cal(C)) = cases(
         -Delta(A, B) / (2 sigma^2) & "if" sigma^2 "is known",
         -(n p) / 2 log(1 + Delta(A, B) / "ESS"(cal(C))) & "if" sigma^2 "is pooled MLE."
       )
     $
     Therefore, Ward's greedy merge is identically the maximum-likelihood agglomerative step for homoscedastic spherical Gaussian mixtures.
]

#proof[
  *Part 1:* Differentiating $f(s) = s(M - s) R^2 / M$ gives $f'(s) = (M - 2s) R^2 / M$ and $f''(s) = -2 R^2 / M < 0$, establishing strict concavity with unique maximum at $s = M / 2$ and minima at endpoints $s in {1, M - 1}$.

  *Part 2:* For a partition $cal(C)$, maximizing the Gaussian log-likelihood over centroids yields the sample means $hat(mu)_C = mu_C$ and residual sum of squares $"ESS"(cal(C))$. When $sigma^2$ is fixed, the profile log-likelihood is $ell(cal(C)) = -(n p / 2) log(2 pi sigma^2) - "ESS"(cal(C)) / (2 sigma^2)$, so merging $A$ and $B$ changes $ell$ by $-("ESS"(cal(C)') - "ESS"(cal(C))) / (2 sigma^2) = -Delta(A, B) / (2 sigma^2)$. When $sigma^2$ is unknown and estimated by its pooled MLE $hat(sigma)^2 = "ESS"(cal(C)) / (n p)$, substituting $hat(sigma)^2$ gives $ell(cal(C)) = -(n p / 2)(log(2 pi "ESS"(cal(C)) / (n p)) + 1)$, so
  $
    ell(cal(C)') - ell(cal(C)) = -(n p) / 2 log( ("ESS"(cal(C)) + Delta(A, B)) / "ESS"(cal(C)) ).
  $
  Because $u |-> log(1 + u)$ is strictly increasing on $[0, infinity)$, maximizing $ell(cal(C)')$ is equivalent to minimizing $Delta(A, B)$.
]

#subsec("D", "Constructive Greedy ESS Suboptimality vs. Optimal k-Means")

Because Ward optimizes the exact same sum-of-squares loss as $k$-means @lloyd1982least (which is NP-hard in general @aloise2009np) subject to hierarchical nesting, a fundamental question is how badly greedy agglomeration can lose to optimal flat $k$-means.

#theorem(title: "Constructive 2-Approximation Suboptimality Even in 1D", origin: "Original Counterexample")[
  For any $epsilon in (0, 1 / 3)$, consider the four-point 1D dataset $X = {0, \, 1, \, 2 - epsilon, \, 3 - epsilon} subset bb(R)$. Then exact Ward clustering yields a $2$-partition with
  $
    "ESS"_"Ward" / "ESS"^(star) = 2 - 2 epsilon + 2 / 3 epsilon^2 -> 2 quad "as" quad epsilon -> 0^(+),
  $
  where $"ESS"^(star) = 1$ is the globally optimal $2$-means cost.
]

#proof[
  At Step 1, all singleton merge costs are $Delta({x_i}, {x_j}) = (1 / 2)(x_i - x_j)^2$. The adjacent costs are $Delta({0}, {1}) = 1 / 2$, $Delta({1}, {2 - epsilon}) = (1 / 2)(1 - epsilon)^2 < 1 / 2$, and $Delta({2 - epsilon}, {3 - epsilon}) = 1 / 2$, while non-adjacent pairs exceed $1$. Hence Ward *must* merge $A = {1, 2 - epsilon}$ at Step 1 with $n_A = 2$, centroid $mu_A = (3 - epsilon) / 2$, and $Delta_1 = (1 / 2)(1 - epsilon)^2$.

  At Step 2, the active clusters are ${0}$, $A$, and ${3 - epsilon}$, with merge costs
  $
    Delta({0}, A) = Delta(A, {3 - epsilon}) &= 2 / 3 ((3 - epsilon) / 2)^2 \
    &= (3 - epsilon)^2 / 6,
  $
  strictly smaller than $Delta({0}, {3 - epsilon}) = (3 - epsilon)^2 / 2$. Thus Step 2 merges either $0$ or $3 - epsilon$ into $A$ with $Delta_2 = (3 - epsilon)^2 / 6$. By Theorem 1 (Part 2), the resulting $2$-clustering achieves
  $
    "ESS"_"Ward" = Delta_1 + Delta_2 &= 1 / 2 (1 - epsilon)^2 + 1 / 6 (3 - epsilon)^2 \
    &= 2 - 2 epsilon + 2 / 3 epsilon^2.
  $
  Meanwhile, the optimal $2$-means partition ${{0, 1}, {2 - epsilon, 3 - epsilon}}$ has centroids $1 / 2$ and $5 / 2 - epsilon$ and cost $"ESS"^(star) = 1 / 2 + 1 / 2 = 1$, proving $"ESS"_"Ward" / "ESS"^(star) -> 2$ as $epsilon -> 0^(+)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational and Memory Profile:* Given a precomputed squared Euclidean distance matrix ($O(n^2 p)$ preprocessing), NN-Chain Ward runs in exact $Theta(n^2)$ time and $Theta(n^2)$ memory @mullner2011modern. Alternatively, in low-to-moderate dimensions $p$, storing only active centroids $mu_C in bb(R)^p$ and sizes $n_C$ evaluates $Delta(A, B)$ on the fly in $O(p)$ time, reducing memory to $Theta(n p)$ while running in $O(n^2 p)$ time. When $n$ is massive, BIRCH @zhang1996birch summarizes points into $m << n$ Clustering Feature (CF) micro-clusters and executes Ward exactly on the CF centroids in $O(n p + m^2)$ time.

*Geometric Inductive Bias:* Ward is invariant to translations and orthogonal rotations of $bb(R)^p$ and equivariant to uniform scaling ($Delta(alpha X) = alpha^2 Delta(X)$), but *not* invariant to anisotropic feature rescaling — feature standardization is essential when units differ. Geometrically, it strongly favors compact, isotropic, equal-variance, balanced-size spherical clusters (Proposition 1) and resists chaining, while its quadratic Euclidean penalty amplifies sensitivity to extreme outliers and causes elongated or non-convex manifolds to be partitioned into compact convex pieces.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$Theta(n^2)$ via NN-Chain on $D^(2)$ (+ $O(n^2 p)$ setup); $O(n^2 p)$ centroid-based.],
    [Space Complexity], [$Theta(n^2)$ matrix mode; $Theta(n p)$ via on-the-fly centroid NN-Chain.],
    [Determinism], [Fully deterministic (up to exact ties); eliminates $k$-means seed variance.],
    [Invariances], [Translation, orthogonal rotation, uniform scaling; *sensitive* to feature scales.],
    [Cluster Geometry], [Compact, isotropic, equal-variance spherical Voronoi-like partitions.],
    [Monotonicity], [Strictly inversion-free by Ward Reducibility (Theorem 2); valid ultrametric.],
    [Outlier Sensitivity], [Immune to chaining; quadratic ESS penalty makes isolated outliers costly.],
    [Primary Bottleneck], [Euclidean requirement (no arbitrary non-metric input) and $O(n^2)$ scaling.],
  ),
  caption: [Operational and theoretical profile of Ward's Minimum-Variance Method.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Deterministic Warm-Start Initialization for $k$-Means and GMMs:* Because Ward directly optimizes the $k$-means ESS objective and the spherical Gaussian log-likelihood (Proposition 1), cutting Ward's dendrogram at level $k$ provides a high-quality, seed-free initialization for Lloyd's $k$-means @lloyd1982least and Expectation-Maximization in `mclust` @fraley2002mclust.
2. *Customer Segmentation and Behavioral Market Profiling:* Commercial analytics datasets typically consist of standardized continuous demographic, recency, frequency, and monetary features where stakeholders require compact, homogeneous, similarly sized segments. Ward's harmonic size-balancing bias and variance minimization match this inductive bias directly.
3. *Functional MRI (fMRI) Brain Parcellation:* Neuroimaging pipelines partition tens of thousands of cortical voxels into functional regions using spatially constrained Ward clustering on BOLD time-series profiles. Enforcing spatial adjacency graphs accelerates Ward to $O(n log n)$ while guaranteeing anatomically contiguous, low-variance cortical parcels.
4. *Geostatistical Regionalization and Econometric Zoning:* Regional planners group census tracts into homogeneous socio-economic districts under contiguity constraints; Ward's exact ANOVA telescoping (Theorem 1) allows analysts to select the number of regions directly from the percentage of total variance explained ($R^2 = 1 - "ESS"(cal(C)_k) / "ESS"({X})$).
5. *Two-Stage Scalable Clustering (BIRCH + Ward):* In high-volume telemetry and industrial IoT pipelines, BIRCH @zhang1996birch compresses millions of points into a compact CF-tree whose sufficient statistics $(n_C, sum x, sum norm(x)_2^2)$ supply Ward's exact merge increment $Delta(A, B)$ with zero information loss relative to exact micro-cluster ESS.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Ward's method is mathematically ill-posed on non-Euclidean dissimilarities (unless embedded via a positive semidefinite kernel) and degrades sharply when datasets contain heavy-tailed outliers (which distort centroids quadratically), severe class-size imbalances (Proposition 1 splits large true clusters prematurely), or intertwined non-convex shapes.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Customer / Market Profiling], [Balanced, low-variance Euclidean segments], [Standardize features; check outliers],
    [fMRI / Spatial Zoning], [Exact $R^2$ variance budget + graph sparsity], [Requires Euclidean / correlation space],
    [GMM / $k$-Means Warm Start], [Matches spherical Gaussian MLE + ESS loss], [Refine cut via Lloyd's $k$-means @lloyd1982least],
    [Non-Euclidean Ecology], [Invalid without Euclidean / RKHS embedding], [Use UPGMA (Alg. 58) or PAM @kaufman1990finding],
    [Curved / Ring Manifolds], [Spherical bias cuts non-convex structures], [Use SLINK (Alg. 56) or HDBSCAN\*],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Ward's Method.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Ward's Method bridges hierarchical clustering, flat $k$-means, kernel methods, and ANOVA internal validity indices through three exact structural links.

#proposition(title: "Calinski--Harabasz Optimality, Kernel Ward, and Hybrid Ward+k-Means", origin: "Original Derivation")[
  1. *Calinski--Harabasz step optimality:* At every coarsening from $k + 1$ to $k$ clusters, Ward's merge maximizes the Calinski--Harabasz variance-ratio index @calinski1974dendrite $"CH"(cal(C)_k) = (B_k / (k - 1)) / ("ESS"(cal(C)_k) / (n - k))$, where $B_k = "ESS"({X}) - "ESS"(cal(C)_k)$ is between-cluster scatter.
  2. *Exact Kernel Ward:* Given any positive semidefinite kernel $kappa(x, y) = chevron.l Phi(x), Phi(y) chevron.r_cal(H)$ on a Hilbert space $cal(H)$, initializing $D_W ({x}, {y}) = kappa(x, x) + kappa(y, y) - 2 kappa(x, y)$ and running Ward's Lance--Williams recurrence (Algorithm Box) executes exact Ward minimum-variance clustering in the implicit feature space $cal(H)$, mirroring Kernel $k$-Means @dhillon2004kernel.
  3. *Monotone hybrid refinement:* Initializing Lloyd's $k$-means @lloyd1982least at Ward's $k$-cluster centroids produces a partition $cal(C)_"hyb"$ satisfying $"ESS"(cal(C)_"hyb") <= "ESS"(cal(C)_"Ward")$, strictly eliminating greedy nesting loss (Theorem 3) whenever a lower-ESS local reassignment exists.
]

#proof[
  *Part 1:* Because $"ESS"({X})$ is constant and $B_k = "ESS"({X}) - "ESS"(cal(C)_k)$, for fixed $k$ the index $"CH"(cal(C)_k) = ((n - k) / (k - 1)) ("ESS"({X}) / "ESS"(cal(C)_k) - 1)$ is strictly decreasing in $"ESS"(cal(C)_k)$. Minimizing $Delta(A, B) = "ESS"(cal(C)_k) - "ESS"(cal(C)_(k+1))$ therefore uniquely maximizes $"CH"(cal(C)_k)$ over all single-merge coarsenings of $cal(C)_(k+1)$.

  *Part 2:* In any real Hilbert space $cal(H)$, the parallel-axis identity and Stewart's centroid identity (Theorem 1) depend solely on inner-product linearity, so $D_W ({x}, {y}) = norm(Phi(x) - Phi(y))_cal(H)^2 = kappa(x, x) + kappa(y, y) - 2 kappa(x, y)$ suffices to run the exact Lance--Williams recurrence without explicit feature coordinates $Phi(x)$.

  *Part 3:* Each Lloyd assignment and centroid-update step non-increases ESS @lloyd1982least, so starting from $cal(C)_"Ward"$ guarantees $"ESS"(cal(C)_"hyb") <= "ESS"(cal(C)_"Ward")$ (and on the Theorem 3 dataset, a single Lloyd iteration immediately recovers the globally optimal partition with $"ESS"^(star) = 1$).
]

Within the Lance--Williams quartet (Algs. 56--59), Ward represents the variance-minimizing pole, contrasting with Single Linkage's connectivity pole, Complete Linkage's diameter pole, and UPGMA's mean-dissimilarity compromise.

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized four-dataset suite, Ward's Method dominates across nearly every benchmark metric: *Iris* $"Sil" = 0.441$, $"ARI" = 0.696$; *Customers* $"Sil" = 0.347$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.292$ ($"ARI" = "N/A"$); with a mean runtime of *2.344 ms*. On *Iris*, Ward achieves by far the highest ground-truth recovery ($"ARI" = 0.696$, outperforming Complete $0.570$ by $+0.126$ and both SLINK and UPGMA $0.560$ by $+0.136$), directly confirming Proposition 1: its homoscedastic Gaussian variance criterion cleanly separates the overlapping _versicolor_ and _virginica_ species even when global Silhouette ($0.441$) slightly favors coarser cuts ($0.489$). On *Customers*, Ward achieves the highest Silhouette ($0.347$, edging out UPGMA $0.334$, Complete $0.323$, and crushing SLINK $0.067$), reflecting its balanced-variance fit for continuous commercial data. On *Multi-Shape*, Ward achieves the highest Silhouette ($0.712$, ahead of SLINK $0.609$, UPGMA $0.452$, and Complete $0.420$), a subtle reminder that Silhouette itself rewards compact low-SSE Euclidean partitions on separated convex sub-blobs even without ground-truth ARI. On *High-D*, Ward again leads all four linkages ($0.292$ vs. $0.245$--$0.257$), demonstrating that centroid-based variance aggregation resists high-dimensional pairwise noise more effectively than extremal or uncentered linkages. Finally, its mean runtime of $2.344$ ms virtually ties SLINK ($2.343$ ms) as fastest — over $2.3 times$ faster than Complete linkage ($5.391$ ms) — confirming the practical speed of exact $Theta(n^2)$ NN-Chain Ward.
