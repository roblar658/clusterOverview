#import "../template/ieee.typ": *

#alg_entry(
  id: 28,
  name: "W-K-Means (Weighted K-Means)",
  category: "Cat. I: Centroid & Partitioning",
  year: "2005",
  paradigm: "Three-Block Alternating Feature-Weighted Minimization",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

In high-dimensional Euclidean datasets, uninformative noise coordinates inflate pairwise distances equally across all clusters, eroding distance contrast @beyer1999nearest and degrading isotropic $k$-means (Alg. 1) @lloyd1982least. *W-$k$-Means* (Automated Variable-Weighting $k$-Means), proposed by Huang, Ng, Rong, and Li @huang2005automated, remedies this weakness by incorporating a continuous simplex vector of feature weights directly into the within-cluster dispersion objective and updating those weights in closed form based on each variable's within-cluster variance.

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ with coordinate entries $x_(i r)$ ($r = 1, dots, d$), let $k in {1, dots, n}$, and let $beta in bb(R) backslash [0, 1]$ be the weight-regularization exponent (standard choice $beta > 1$). Denote the probability simplex by $Delta_(d - 1) := {bold(w) in bb(R)_(>= 0)^d : sum_(r=1)^d w_r = 1}$. W-$k$-Means seeks a disjoint $k$-partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$, centroids $bold(mu) = (mu_1, dots, mu_k) in (bb(R)^d)^k$, and feature weights $bold(w) in Delta_(d - 1)$ minimizing the *weighted within-cluster sum of squares* @huang2005automated:
$ J_beta (cal(C), bold(mu), bold(w)) &= sum_(j=1)^k sum_(i in C_j) sum_(r=1)^d w_r^beta (x_(i r) - mu_(j r))^2 \
  &= sum_(r=1)^d w_r^beta D_r (cal(C), bold(mu)), $
where $D_r (cal(C), bold(mu)) := sum_(j=1)^k sum_(i in C_j) (x_(i r) - mu_(j r))^2 >= 0$ is the total within-cluster dispersion along coordinate $r$. Intuitively, for $beta > 1$, coordinates that exhibit tight within-cluster variance (small $D_r$) receive larger weight $w_r$, whereas diffuse noise features (large $D_r$) are smoothly attenuated.

#definition(title: "Feature-Weighted Centroidal Equilibrium")[
  A triple $(cal(C), bold(mu), bold(w))$ with $bold(w) in Delta_(d - 1)$ is a *feature-weighted centroidal equilibrium* of $J_beta$ if: (i) every $i in C_j$ satisfies $sum_(r=1)^d w_r^beta (x_(i r) - mu_(j r))^2 <= sum_(r=1)^d w_r^beta (x_(i r) - mu_(ell r))^2$ for all $ell$ (weighted Voronoi optimality); (ii) $mu_j = |C_j|^(-1) sum_(i in C_j) x_i$ for each non-empty $C_j$ (centroid optimality); and (iii) $bold(w) in arg min_(u in Delta_(d-1)) sum_(r=1)^d u_r^beta D_r (cal(C), bold(mu))$ (simplex weight optimality).
]

Huang et al. @huang2005automated solve this non-convex three-block program by alternating exact minimizations over $cal(C)$, $bold(mu)$, and $bold(w)$ (Algorithm 1), generalizing both Euclidean $k$-means and categorical $k$-modes @huang1998extensions within a unified variable-weighting framework.

#algbox(
  title: "W-K-Means Three-Block Alternating Minimization",
  inputs: [Standardized dataset $cal(X) subset.eq bb(R)^d$, clusters $k >= 1$, exponent $beta > 1$, initial centroids $bold(mu)^((0))$, initial weights $bold(w)^((0)) = (1 / d, dots, 1 / d)$, max iterations $T$.],
  outputs: [Partition $cal(C)$, centroids $bold(mu)$, and learned feature weights $bold(w) in Delta_(d-1)$.],
)[
  1. *For* $t = 0, 1, dots, T - 1$ *do*:
  2. #h(1em) *Weighted Voronoi assignment:* For each $i = 1, dots, n$, set
     $ c_i^((t)) = arg min_(1 <= j <= k) sum_(r=1)^d (w_r^((t)))^beta (x_(i r) - mu_(j r)^((t)))^2, $
     breaking ties by smallest cluster index $j$; form $C_j^((t)) = {i : c_i^((t)) = j}$.
  3. #h(1em) *Centroid update:* For each non-empty $C_j^((t))$, compute
     $ mu_j^((t+1)) = 1 / (|C_j^((t))|) sum_(i in C_j^((t))) x_i. $
  4. #h(1em) *Coordinate dispersions:* For each $r = 1, dots, d$, compute
     $ D_r^((t+1)) = sum_(j=1)^k sum_(i in C_j^((t))) (x_(i r) - mu_(j r)^((t+1)))^2. $
  5. #h(1em) *Feature-weight update:* Let $Z = {r : D_r^((t+1)) = 0}$. If $Z != emptyset$, set $w_r^((t+1)) = |Z|^(-1)$ for $r in Z$ and $0$ otherwise; else set
     $ w_r^((t+1)) = ( (D_r^((t+1)))^(-1 / (beta - 1)) ) / ( sum_(q=1)^d (D_q^((t+1)))^(-1 / (beta - 1)) ). $
  6. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$, *break*.
  7. *Return* $(cal(C)^((t)), bold(mu)^((t+1)), bold(w)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five core theoretical properties of W-$k$-Means: (A) the exact closed-form KKT solution and phase transition across all regimes of $beta$; (B) weight-invariance of centroid projection and strict three-block finite-step convergence; (C) an explicit counterexample proving orthogonal non-invariance and catastrophic unstandardized scale collapse; and (D) an analytical high-dimensional signal-to-noise weight bound dictating the admissible range of $beta$.

#subsec("A", "Closed-Form KKT Weights and Phase Transition in beta")

#theorem(title: "Exact Simplex Weight Optimality and Phase Regimes in beta", origin: "Original Proof")[
  Fix $(cal(C), bold(mu))$ with coordinate dispersions $D_r := D_r (cal(C), bold(mu)) >= 0$ and zero-dispersion set $Z = {r : D_r = 0}$. Consider $min_(bold(w) in Delta_(d-1)) f_beta (bold(w))$ where $f_beta (bold(w)) := sum_(r=1)^d w_r^beta D_r$:
  1. *Strictly convex informative regime ($beta > 1$):* If $Z = emptyset$, the unique global minimizer is
     $ w_r^* = ( D_r^(-1 / (beta - 1)) ) / ( sum_(q=1)^d D_q^(-1 / (beta - 1)) ), $
     achieving minimum value $f_beta (bold(w)^*) = (sum_(r=1)^d D_r^(-1 / (beta - 1)))^(1 - beta)$. If $Z != emptyset$, $min f_beta = 0$, attained by every $bold(w) in Delta_(d - 1)$ supported on $Z$.
  2. *Concave degenerate regime ($0 < beta <= 1$):* Every global minimizer lies on a simplex vertex $e_(r^*)$ with $r^* in arg min_r D_r$, collapsing the metric onto a single coordinate.
  3. *Inverted-weighting regime ($beta < 0$, $Z = emptyset$):* The unique minimizer is $w_r^* = D_r^(1 / (1 - beta)) / sum_(q=1)^d D_q^(1 / (1 - beta))$, which assigns *strictly larger* weight to noisier (higher-$D_r$) coordinates.
]

#proof[
  For any coordinate with $D_r > 0$, the second derivative of $w_r |-> w_r^beta D_r$ on $(0, 1]$ is $beta (beta - 1) w_r^(beta - 2) D_r$.

  (1) When $beta > 1$ and $Z = emptyset$, $beta (beta - 1) > 0$, so $f_beta$ is strictly convex on $Delta_(d - 1)$. Introducing Lagrange multiplier $lambda in bb(R)$ for $sum_r w_r = 1$ and dual multipliers $nu_r >= 0$ for $w_r >= 0$, the KKT conditions are
  $ beta w_r^(beta - 1) D_r = lambda + nu_r, quad nu_r w_r = 0, quad r = 1, dots, d. $
  At least one coordinate has $w_(r_0) > 0$, forcing $nu_(r_0) = 0$ and hence $lambda = beta w_(r_0)^(beta - 1) D_(r_0) > 0$. If any coordinate had $w_r = 0$, the left-hand side would vanish while $lambda + nu_r >= lambda > 0$, a contradiction. Thus $w_r > 0$ and $nu_r = 0$ for every $r$, giving $w_r = (lambda / (beta D_r))^(1 / (beta - 1))$. Enforcing $sum_(r=1)^d w_r = 1$ yields $w_r^* = S^(-1) D_r^(-1 / (beta - 1))$ where $S := sum_(q=1)^d D_q^(-1 / (beta - 1))$. Substituting back gives
  $ f_beta (bold(w)^*) &= sum_(r=1)^d S^(-beta) D_r^(-beta / (beta - 1)) D_r \
    &= S^(-beta) sum_(r=1)^d D_r^(-1 / (beta - 1)) = S^(1 - beta). $
  If $Z != emptyset$, $f_beta (bold(w)) >= 0$ everywhere and vanishes if and only if $w_r = 0$ for all $r in.not Z$.

  (2) For $0 < beta <= 1$, $beta (beta - 1) <= 0$, so $f_beta$ is concave on the compact simplex $Delta_(d - 1)$. The minimum of a concave function over a polytope is attained at an extreme point $e_r$, where $f_beta (e_r) = D_r$, hence at $r^* in arg min_r D_r$.

  (3) For $beta < 0$ and $Z = emptyset$, $beta (beta - 1) > 0$, so $f_beta$ is strictly convex on the open simplex with barrier $w_r^beta D_r -> +infinity$ as $w_r -> 0^+$. Stationarity $beta w_r^(beta - 1) D_r = lambda < 0$ yields $w_r^(1 - beta) = (-beta / lambda) D_r$, so $w_r^* prop D_r^(1 / (1 - beta))$ is strictly increasing in $D_r$, contradicting noise suppression.
]

#subsec("B", "Three-Block Lyapunov Descent and Finite Termination")

#lemma(title: "Weight-Invariance of Optimal Centroid Projection", origin: "Original Proof")[
  For any fixed partition $cal(C)$ and any weights $bold(w) in Delta_(d - 1)$ with $beta > 1$, the unweighted arithmetic centroid $mu_j^* = |C_j|^(-1) sum_(i in C_j) x_i$ globally minimizes $mu_j |-> sum_(i in C_j) sum_(r=1)^d w_r^beta (x_(i r) - mu_(j r))^2$, with exact reduction
  $ sum_(i in C_j) sum_(r=1)^d w_r^beta (x_(i r) - mu_(j r))^2 &= sum_(i in C_j) sum_(r=1)^d w_r^beta (x_(i r) - mu_(j r)^*)^2 \
    &quad + |C_j| sum_(r=1)^d w_r^beta (mu_(j r) - mu_(j r)^*)^2. $
]

#proof[
  For each coordinate $r in {1, dots, d}$, the 1D parallel-axis identity gives $sum_(i in C_j) (x_(i r) - mu_(j r))^2 = sum_(i in C_j) (x_(i r) - mu_(j r)^*)^2 + |C_j| (mu_(j r) - mu_(j r)^*)^2$, because $sum_(i in C_j) (x_(i r) - mu_(j r)^*) = 0$. Multiplying by $w_r^beta >= 0$ and summing over $r = 1, dots, d$ establishes the identity directly; non-negativity of the second sum confirms global optimality of $mu_j^*$ independently of $bold(w)$.
]

#theorem(title: "Strict Three-Block Descent and Finite-Step Convergence", origin: "Classical — Cited")[
  For any $beta > 1$, deterministic smallest-index Voronoi tie-breaking, and non-empty clusters, Algorithm 1 generates a sequence ${(cal(C)^((t)), bold(mu)^((t)), bold(w)^((t)))}$ satisfying @huang2005automated:
  1. $J_beta (cal(C)^((t+1)), bold(mu)^((t+1)), bold(w)^((t+1))) <= J_beta (cal(C)^((t)), bold(mu)^((t)), bold(w)^((t)))$ at every step, with exact centroid drop given by Lemma 1.
  2. Algorithm 1 terminates in at most $T <= sum_(m=1)^k S(n, m) <= k^n / k!$ iterations at a feature-weighted centroidal equilibrium.
]

#proof[
  At iteration $t$, Step 2 minimizes $J_beta (cal(C), bold(mu)^((t)), bold(w)^((t)))$ over $cal(C)$ pointwise for each $x_i$; Step 3 minimizes $J_beta (cal(C)^((t)), bold(mu), bold(w)^((t)))$ over $bold(mu)$ by Lemma 1; and Step 5 minimizes $J_beta (cal(C)^((t)), bold(mu)^((t+1)), bold(w))$ over $bold(w) in Delta_(d - 1)$ by Theorem 1(1). Chaining the three exact block minimizations establishes non-increasing monotonicity of $J_beta$.

  Crucially, despite the continuous weight vector $bold(w) in Delta_(d - 1)$, both $bold(mu)^((t+1))$ and $bold(w)^((t+1))$ are *uniquely determined* by the discrete partition $cal(C)^((t))$ alone: $bold(mu)^((t+1)) = bold(mu)^* (cal(C)^((t)))$ via arithmetic means, and $bold(w)^((t+1)) = bold(w)^* (cal(C)^((t)))$ via Step 5 (uniquely by strict convexity when $Z = emptyset$, and deterministically uniform over $Z$ when $Z != emptyset$). Define the reduced combinatorial energy $J_beta^* (cal(C)) := J_beta (cal(C), bold(mu)^* (cal(C)), bold(w)^* (cal(C)))$. Suppose $J_beta^* (cal(C)^((t))) = J_beta^* (cal(C)^((t-1)))$ at some iteration $t >= 1$. Then all three substeps incur zero decrease:
  - If $J_beta^* (cal(C)^((t-1))) > 0$, then $Z = emptyset$ and $w_r^((t)) > 0$ for all $r$. Zero decrease in Step 3 (Lemma 1) forces $bold(mu)^((t+1)) = bold(mu)^((t))$, and strict convexity of Step 5 (Theorem 1(1)) forces $bold(w)^((t+1)) = bold(w)^((t))$. Consequently, Step 2 at iteration $t + 1$ evaluates the identical parameter pair $(bold(mu)^((t)), bold(w)^((t)))$, and deterministic smallest-index tie-breaking reproduces $cal(C)^((t+1)) = cal(C)^((t))$, halting immediately.
  - If $J_beta^* (cal(C)^((t-1))) = 0$, then $Z != emptyset$ and every point already has zero weighted distance on $Z$. Step 2 only reassigns ties to the smallest cluster index matching each point's $Z$-coordinates, leaving $mu_(j r)$ ($r in Z$), $Z$, and $bold(w)$ invariant, so iteration $t + 1$ reproduces $cal(C)^((t+1)) = cal(C)^((t))$ and halts.

  Hence every iteration prior to the terminal tie-resolution check strictly decreases $J_beta^* (cal(C))$, so no partition can recur. Because the number of partitions of $n$ points into at most $k$ non-empty clusters is $sum_(m=1)^k S(n, m) <= k^n / k!$, Algorithm 1 halts in finitely many steps at a feature-weighted centroidal equilibrium.
]

#subsec("C", "Orthogonal Non-Invariance and Scale-Collapse Counterexample")

Whereas Lloyd's $k$-means is invariant under arbitrary rotations $Q in O(d)$, axis-aligned variable weighting breaks rotational symmetry and amplifies sensitivity to feature scaling.

#proposition(title: "Rotational Non-Invariance and Unstandardized Scale Collapse", origin: "Original Counterexample")[
  Let $beta > 1$. Then:
  1. *Rotational non-invariance:* W-$k$-Means is invariant under signed coordinate permutations and translations, but *not* under general orthogonal rotations $Q in O(d)$.
  2. *Catastrophic scale collapse:* Consider a 2D dataset where coordinate $1$ has informative cluster structure with minimal attainable dispersion $D_(1, "min") > 0$, while coordinate $2$ is pure uninformative noise with total variance $sigma_2^2 > 0$ rescaled by $gamma > 0$. For *every* partition $cal(C)$, as $gamma -> 0^+$, the learned weights collapse onto the noise coordinate:
     $ w_1^* (cal(C)) / (w_2^* (cal(C))) <= ( (n gamma^2 sigma_2^2) / D_(1, "min") )^(1 / (beta - 1)) -> 0, quad w_2^* -> 1. $
]

#proof[
  (1) Construct 4 points in $bb(R)^2$: $(0, -1), (0, 1), (10, -1), (10, 1)$ with $k = 2$ and $beta = 2$. Along the canonical axes, the vertical-pair partition yields $D_1 = 0, D_2 = 4$, giving $bold(w)^* = (1, 0)$ and $J_2 = 0$. Rotating all points by $45 degree$ via $Q = 2^(-1/2) mat(1, -1; 1, 1)$ maps each vertical pair's difference vector $(0, 2)^T$ to $(-sqrt(2), sqrt(2))^T$, making $D_1 = D_2 = 2 > 0$, which forces uniform weights $bold(w)^* = (1 / 2, 1 / 2)$ and $J_2 = 1 > 0$. Thus orthogonal rotations destroy axis-aligned sparsity.

  (2) Under any partition $cal(C)$, within-cluster dispersion never exceeds total sample dispersion, so $D_2 (cal(C)) <= n gamma^2 sigma_2^2$, whereas $D_1 (cal(C)) >= D_(1, "min") > 0$. By Theorem 1(1),
  $ w_1^* / w_2^* &= ( D_2 (cal(C)) / (D_1 (cal(C))) )^(1 / (beta - 1)) \
    &<= ( (n gamma^2 sigma_2^2) / D_(1, "min") )^(1 / (beta - 1)) = O(gamma^(2 / (beta - 1))). $
  As $gamma -> 0^+$, $w_1^* / w_2^* -> 0$, forcing $w_2^* -> 1$ and $w_1^* -> 0$—assigning nearly $100\%$ weight to pure noise merely because its units are numerically smaller.
]

#subsec("D", "High-Dimensional Signal-to-Noise Weight Ratio Bound")

Standardizing every feature to unit variance eliminates the scale trap of Proposition 1(2) and allows an exact quantitative analysis of noise-feature suppression.

#proposition(title: "Signal-to-Noise Subspace Concentration Bound", origin: "Original Derivation")[
  Suppose every feature $r in {1, dots, d}$ is standardized so that $sum_(i=1)^n (x_(i r) - overline(x)_r)^2 = n$. Decompose ${1, dots, d} = cal(S) union cal(N)$ into $d_S = |cal(S)|$ informative coordinates where the true partition $cal(C)^*$ explains variance fraction $R_r^2 >= R_"min"^2 > 0$ (i.e., $D_r (cal(C)^*) = n (1 - R_r^2)$) and $d_N = d - d_S$ noise coordinates satisfying $D_q (cal(C)^*) >= n (1 - delta)$ with $0 <= delta < R_"min"^2$. Then at $cal(C)^*$:
  1. Every signal-to-noise weight ratio satisfies $w_r^* / w_q^* >= ((1 - delta) / (1 - R_r^2))^(1 / (beta - 1))$.
  2. Total weight captured by the signal subspace $cal(S)$ satisfies
     $ sum_(r in cal(S)) w_r^* >= [ 1 + d_N / d_S ( (1 - R_"min"^2) / (1 - delta) )^(1 / (beta - 1)) ]^(-1). $
  3. Consequently, ensuring $sum_(r in cal(S)) w_r^* >= 1 / 2$ when $d_N > d_S$ requires
     $ 1 < beta <= 1 + ( log ((1 - delta) / (1 - R_"min"^2)) ) / ( log (d_N / d_S) ). $
]

#proof[
  For every $r in cal(S)$, $D_r (cal(C)^*) <= n (1 - R_"min"^2)$, while for every $q in cal(N)$, $D_q (cal(C)^*) >= n (1 - delta)$. Because $beta > 1$, the map $u |-> u^(-1 / (beta - 1))$ is strictly decreasing, giving
  $ w_r^* / w_q^* &= ( D_q (cal(C)^*) / (D_r (cal(C)^*)) )^(1 / (beta - 1)) \
    &>= ( (1 - delta) / (1 - R_r^2) )^(1 / (beta - 1)), $
  proving Claim 1. Summing the numerator and denominator separately over $cal(N)$ and $cal(S)$ bounds the noise-to-signal weight ratio:
  $ ( sum_(q in cal(N)) w_q^* ) / ( sum_(r in cal(S)) w_r^* ) &= ( sum_(q in cal(N)) D_q^(-1 / (beta - 1)) ) / ( sum_(r in cal(S)) D_r^(-1 / (beta - 1)) ) \
    &<= d_N / d_S ( (1 - R_"min"^2) / (1 - delta) )^(1 / (beta - 1)). $
  Using $sum_(r in cal(S)) w_r^* + sum_(q in cal(N)) w_q^* = 1$ immediately yields Claim 2. Setting the lower bound in Claim 2 $>= 1 / 2$ is equivalent to $(d_N / d_S) ((1 - R_"min"^2) / (1 - delta))^(1 / (beta - 1)) <= 1$; taking logarithms and solving for $beta - 1 > 0$ proves Claim 3.
]

#remark(title: "Breakdown Point and NP-Hardness")[
  Because centroids in Step 3 remain unweighted sample means, corrupting a single point $x_1 -> infinity$ diverges both its centroid and $D_r$, so W-$k$-Means inherits breakdown point $epsilon^* = 1 / n$ (and NP-hardness @aloise2009np from the fixed-$bold(w)$ subproblem).
]

#sec("III", "Complexity, Geometry and Operational Profile")

W-$k$-Means retains the exact asymptotic complexity of Lloyd's algorithm: by pre-scaling feature columns by $w_r^(beta / 2)$ in $O((n + k) d)$ time, Step 2 reduces to a standard Euclidean BLAS-3 distance matrix product costing $Theta(n k d)$. Step 3 centroid accumulation and Step 4 coordinate dispersion accumulation each take $O(n d)$ time, and Step 5 simplex normalization takes $O(d)$ time. Overall complexity is $O(n k d T)$ with $O((n + k) d)$ memory.

Geometrically, W-$k$-Means induces an *axis-aligned Mahalanobis Voronoi tessellation* governed by the diagonal metric tensor $W_beta = "diag"(w_1^beta, dots, w_d^beta)$. Bisector hyperplanes $(x - (mu_j + mu_ell) / 2)^T W_beta (mu_j - mu_ell) = 0$ tilt away from Euclidean perpendicular bisectors to prioritize low-dispersion separating axes.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for W-K-Means],
    [Time / Space Complexity], [$O(n k d T)$ time (same order as Alg. 1 via column scaling); $O((n + k) d)$ memory.],
    [Geometric Inductive Bias], [Axis-aligned ellipsoidal Voronoi cells with shared global diagonal metric $W_beta = "diag"(w_r^beta)$.],
    [Determinism / Seeding], [Deterministic given $(bold(mu)^((0)), bold(w)^((0)))$; sensitive to both centroid seeding and $beta$.],
    [Invariances], [Invariant under signed axis permutations & translations; *not* invariant to rotations or unstandardized scales (Prop. 1).],
    [Outlier Breakdown Point], [$epsilon^* = 1 / n$ (unweighted centroids; also vulnerable to $D_r = 0$ constant-column collapse).],
    [Recommended Upgrades], [Minkowski Weighted (Alg. 29) for cluster-specific $L_p$ weights; PROCLUS (Alg. 136) or COSA (Alg. 145) for per-cluster subspaces.],
  ),
  caption: [Summary of theoretical guarantees and operational profile for W-K-Means (Alg. 28).],
)

#sec("IV", "Engineering and Scientific Applications")

W-$k$-Means is well suited to tabular and scientific datasets that have been Z-score standardized and contain a *global* subset of cluster-informative features mixed with irrelevant nuisance attributes @aggarwal2013data @hastie2009elements.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Transcriptomic & Proteomic Cohort Discovery:* High-throughput biological assays measure hundreds to thousands of genes/proteins, of which only a compact global pathway separates disease subtypes while housekeeping genes fluctuate independently. Soft inverse-dispersion weighting amplifies discriminative biomarkers without brittle hard thresholding, while learned weights $bold(w)$ serve directly as interpretable feature-importance rankings.
2. *Enterprise Customer & Financial Behavioral Segmentation:* Standardized CRM tables combine informative behavioral KPIs (recency, spend elasticity, tenure) with low-signal demographic or operational fields. W-$k$-Means automatically suppresses uninformative columns and preserves interpretable unweighted arithmetic centroids (Lemma 1) for business reporting.
3. *Industrial Condition Monitoring & Multi-Sensor Diagnostics:* Rotating machinery and semiconductor fabs log dozens of heterogeneous telemetry channels (vibrations, temperatures, currents) after Z-score calibration. Fault regimes manifest sharply on a subset of physical transducers; W-$k$-Means concentrates cluster boundaries onto those low-residual sensors while down-weighting ambient drift channels.
4. *Multispectral Remote Sensing & Land-Cover Mapping:* Satellite imagery records reflectance across visible, near-infrared, and thermal bands, some degraded by atmospheric scattering. Adaptive band weighting attenuates noisy spectral channels and sharpens land-cover Voronoi separation at linear $O(n k d)$ throughput.
5. *Structured Text & Tabular Feature Auditing:* On standardized topic-embedding or document-attribute tables, W-$k$-Means simultaneously clusters documents and identifies which latent dimensions exhibit consistent within-topic cohesion versus high intra-topic dispersion.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes govern practical deployment: (a) *missing feature standardization*, which causes immediate weight collapse onto whichever noise feature has the smallest numeric units (Proposition 1(2)); (b) *cluster-specific or rotated subspaces*, because W-$k$-Means learns a single *global* axis-aligned vector $bold(w)$ shared by all $k$ clusters (failing when cluster 1 lives in coordinates ${1, 2}$ and cluster 2 lives in ${3, 4}$—resolved by PROCLUS Alg. 136 @aggarwal1999proclus, ORCLUS Alg. 138 @aggarwal2000orclus, COSA Alg. 145 @friedman2004cosa, or survey methods @kriegel2009subspacesurvey); (c) *near-constant or discrete zero-dispersion columns* ($D_r -> 0$), which hijack $100\%$ of simplex weight unless a small floor $D_r + epsilon$ is added; and (d) *poor tuning of $beta$*, since $beta <= 1$ collapses to 1D (Theorem 1(2)) while excessively large $beta$ violates Proposition 2(3) and fails to suppress high-dimensional noise.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why W-K-Means Fits], [Caution / Alternative],
    [Bioinformatics], [Global biomarker ranking + soft noise suppression], [Z-score first; add $epsilon$ floor],
    [CRM / Tabular], [Suppresses noisy KPIs; keeps mean profiles], [Use Alg. 8 for mixed categories],
    [Sensor Diagnostics], [Locks onto fault-sensitive channels at $O(n k d)$], [Use Alg. 26 if outlier spikes],
    [Multi-Subspace Data], [Fast baseline when global subspace exists], [Use Alg. 29/136/145 per-cluster],
    [Correlated / Oblique], [Axis-aligned Mahalanobis approximation], [Use Alg. 81 (GMM) or Alg. 138],
  ),
  caption: [Domain suitability matrix and engineering guidance for W-K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Eliminating the simplex weight vector $bold(w)$ reveals a striking exact duality linking W-$k$-Means (Alg. 28) simultaneously to Power $k$-Means (Alg. 27) @xu2019power, diagonal Gaussian Mixture Models (Alg. 81 / MCLUST Alg. 93) @fraley2002mclust, and Minkowski Weighted $k$-Means (Alg. 29) @deamorim2012minkowski.

#proposition(title: "Feature-Space Power-Mean Duality and Diagonal-GMM Limit", origin: "Original Unified Proof")[
  For any partition $cal(C)$ and centroids $bold(mu)$ with positive dispersions $D_1, dots, D_d > 0$, let $J_beta^* (cal(C), bold(mu)) := min_(bold(w) in Delta_(d-1)) J_beta (cal(C), bold(mu), bold(w))$ for $beta > 1$, and let $M_s (D_1, dots, D_d) = (d^(-1) sum_(r=1)^d D_r^s)^(1 / s)$ denote the power mean of order $s$. Then:
  1. *Exact Power-Mean Dual:* Setting $s = -1 / (beta - 1) < 0$,
     $ J_beta^* (cal(C), bold(mu)) = d^(1 - beta) M_(-1 / (beta - 1)) (D_1, dots, D_d). $
     In particular, at $beta = 2$, $d dot J_2^*$ equals the exact *harmonic mean* of coordinate dispersions.
  2. *Asymptotic Limits:* As $beta -> 1^+$ ($s -> -infinity$), $d^(beta - 1) J_beta^* -> min_(1 <= r <= d) D_r$. Conversely, as $beta -> +infinity$ ($s -> 0^-$), the rescaled objective converges to the geometric-mean / determinant criterion of a shared-diagonal-covariance Gaussian Mixture Model:
     $ lim_(beta -> +infinity) d^(beta - 1) J_beta^* (cal(C), bold(mu)) = ( product_(r=1)^d D_r )^(1 / d), $
     while the effective metric weight ratios converge to inverse-variance Mahalanobis scaling: $(w_r^* / w_q^*)^beta -> D_q / D_r$.
]

#proof[
  (1) By Theorem 1(1), $J_beta^* (cal(C), bold(mu)) = (sum_(r=1)^d D_r^(-1 / (beta - 1)))^(1 - beta)$. Setting $s = -1 / (beta - 1) < 0$ gives $1 / s = 1 - beta$, so
  $ J_beta^* (cal(C), bold(mu)) &= d^(1 / s) ( 1 / d sum_(r=1)^d D_r^s )^(1 / s) \
    &= d^(1 - beta) M_s (D_1, dots, D_d). $
  At $beta = 2$, $s = -1$ and $d^(beta - 1) = d$, yielding $d dot J_2^* = (d^(-1) sum_(r=1)^d D_r^(-1))^(-1)$.

  (2) As $beta -> 1^+$, $s = -1 / (beta - 1) -> -infinity$, so $M_s (D_1, dots, D_d) -> min_r D_r$. As $beta -> +infinity$, $s -> 0^-$, and the classical zero-order limit of the power mean gives $lim_(s -> 0^-) M_s (D_1, dots, D_d) = (product_(r=1)^d D_r)^(1 / d)$. Taking logarithms yields $d^(-1) sum_(r=1)^d log D_r$, matching the profile negative log-likelihood of a hard-assignment Gaussian mixture with shared diagonal covariance $"diag"(D_1 / n, dots, D_d / n)$ @fraley2002mclust. Finally, from Theorem 1(1),
  $ ( w_r^* / w_q^* )^beta = ( D_q / D_r )^(beta / (beta - 1)) -> D_q / D_r quad (beta -> +infinity), $
  confirming inverse-variance diagonal Mahalanobis weighting in the large-$beta$ limit.
]

Proposition 3 exposes a clean structural symmetry within Category I: Power $k$-Means (Alg. 27) applies a negative power mean across *clusters* $j in {1, dots, k}$ to smooth local traps, whereas W-$k$-Means (Alg. 28) applies a negative power mean across *features* $r in {1, dots, d}$ to learn global feature relevance. Minkowski Weighted $k$-Means (Alg. 29) @deamorim2012minkowski directly extends Alg. 28 by replacing global weights $w_r$ with cluster-specific weights $w_(j r)$ and unifying the distance exponent with $beta$ via $sum_r w_(j r)^beta |x_(i r) - mu_(j r)|^beta$.

In real empirical benchmarks on the ClusterOverview 4-dataset suite @rousseeuw1987silhouettes @hubert1985comparing, W-$k$-Means attains Iris Sil $= 0.478$ (ARI $= 0.441$), Customers Sil $= 0.328$, Multi-Shape Sil $= 0.718$ (ARI N/A), and High-D Subspace+Noise Sil $= 0.243$ (ARI N/A), with a fast mean runtime of $0.456$ ms. Notably, adaptive feature weighting yields the highest Silhouette among Algs. 26–28 on both Multi-Shape ($0.718$) and Iris ($0.478$) by amplifying dominant low-variance coordinates, yet lowers Iris ARI ($0.441$) as down-weighting sepal dimensions shifts overlapping class boundaries. Conversely, its modest score on High-D Subspace+Noise ($0.243$) concretely validates Proposition 2 and Section IV.B: a single *global* weight vector cannot isolate disjoint cluster-specific subspaces when $d_N / d_S$ is large, motivating cluster-specific subspace weighting (Algs. 29, 136, 145).
