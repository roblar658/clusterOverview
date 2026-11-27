#import "../template/ieee.typ": *

#alg_entry(
  id: 1,
  name: "K-Means (Lloyd's Algorithm)",
  category: "Cat. I: Centroid & Partitioning",
  year: "1957 / 1982",
  paradigm: "Alternating Bregman / Voronoi Minimization",
  complexity: $O(n k d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Variational Objective")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset equipped with the Euclidean norm $norm(dot)_2$, and let $k in {1, dots, n}$ denote the target number of clusters. The *Euclidean $k$-means clustering problem* seeks a disjoint partition $cal(C) = {C_1, dots, C_k}$ of ${1, dots, n}$ together with prototype vectors $bold(mu) = (mu_1, dots, mu_k) in (bb(R)^d)^k$ minimizing the *within-cluster sum of squares (WCSS)* functional @lloyd1982least:
$ J(cal(C), bold(mu)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - mu_j)_2^2. $
Minimizing $J$ over $bold(mu)$ for a fixed partition yields the purely combinatorial formulation over partitions alone:
$ J^*(cal(C)) = sum_(j=1)^k sum_(i in C_j) norm(x_i - overline(x)_(C_j))_2^2 = sum_(j=1)^k 1 / (2 |C_j|) sum_(i, ell in C_j) norm(x_i - x_ell)_2^2, $
where $overline(x)_(C_j) = |C_j|^(-1) sum_(i in C_j) x_i$ is the arithmetic centroid of $C_j$. Conversely, fixing $bold(mu)$ induces the *Voronoi tessellation* of $bb(R)^d$ whose cells are closed convex polyhedra:
$ V_j(bold(mu)) = { x in bb(R)^d : norm(x - mu_j)_2^2 <= norm(x - mu_ell)_2^2 space forall ell != j }. $

#definition(title: "Centroidal Voronoi Equilibrium")[
  A pair $(cal(C), bold(mu))$ is a *discrete centroidal Voronoi tessellation (CVT)*—equivalently, a coordinate-wise stationary point of $J$—if (i) every $i in C_j$ satisfies $x_i in V_j(bold(mu))$ (Voronoi optimality), and (ii) for every non-empty cluster $C_j$, $mu_j = |C_j|^(-1) sum_(i in C_j) x_i$ (centroid optimality).
]

Lloyd's algorithm @lloyd1982least solves the non-convex joint optimization via alternating block coordinate descent (Algorithm 1), refining assignments and centroids until a discrete CVT is reached.

#algbox(
  title: "Lloyd's K-Means Algorithm",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, integer $k >= 1$, initial centroids $bold(mu)^((0))$, max iterations $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$ and centroids $bold(mu)$.],
)[
  1. *For* $t = 0, 1, dots, T - 1$ *do*:
  2. #h(1em) *Assignment step (Voronoi projection):* For each $i in {1, dots, n}$, assign
     $ c_i^((t)) = arg min_(j in {1, dots, k}) norm(x_i - mu_j^((t)))_2^2, $
     breaking ties consistently by smallest cluster index $j$.
  3. #h(1em) Set $C_j^((t)) = {i : c_i^((t)) = j}$ for $j = 1, dots, k$.
  4. #h(1em) *Update step (Mean projection):* For each non-empty $C_j^((t))$, compute
     $ mu_j^((t+1)) = 1 / (|C_j^((t))|) sum_(i in C_j^((t))) x_i. $
     (Empty clusters are reseeded to the point maximizing $min_ell norm(x_i - mu_ell^((t)))_2^2$.)
  5. #h(1em) *If* $cal(C)^((t)) = cal(C)^((t-1))$ (or $bold(mu)^((t+1)) = bold(mu)^((t))$), *break*.
  6. *Return* $(cal(C)^((t)), bold(mu)^((t+1)))$.
]

#sec("II", "Formal Mathematical Proofs and Guarantees")

We establish five core theoretical properties governing Lloyd's algorithm: (1) strict energy descent and finite-step termination with a quantitative step-decrease bound; (2) equivalence to an exact block Gauss–Newton step with local quadratic convergence; (3) global NP-hardness, worst-case exponential iteration lower bounds, and an explicit adversarial counterexample proving unbounded suboptimality; (4) geometric invariance and finite-sample breakdown point zero; and (5) strong statistical consistency of empirical $k$-means minimizers.

#subsec("A", "Quantitative Energy Descent and Finite Termination")

#lemma(title: "Parallel-Axis / Bregman Pythagorean Identity", origin: "Original Proof")[
  For any finite non-empty subset $S subset.eq bb(R)^d$ with centroid $overline(x)_S = |S|^(-1) sum_(x in S) x$ and any vector $c in bb(R)^d$,
  $ sum_(x in S) norm(x - c)_2^2 = sum_(x in S) norm(x - overline(x)_S)_2^2 + |S| norm(overline(x)_S - c)_2^2. $
  Consequently, $overline(x)_S$ is the *unique* global minimizer of $c |-> sum_(x in S) norm(x - c)_2^2$.
]

#proof[
  Expand the squared Euclidean norm by inserting $plus.minus overline(x)_S$:
  $ &sum_(x in S) norm((x - overline(x)_S) + (overline(x)_S - c))_2^2 \
    &= sum_(x in S) norm(x - overline(x)_S)_2^2 + |S| norm(overline(x)_S - c)_2^2 \
    &quad + 2 chevron.l sum_(x in S) (x - overline(x)_S), space overline(x)_S - c chevron.r. $
  By definition of the arithmetic mean, $sum_(x in S) (x - overline(x)_S) = sum_(x in S) x - |S| overline(x)_S = bold(0)$. Thus the cross-term vanishes identically, leaving $sum_(x in S) norm(x - overline(x)_S)_2^2 + |S| norm(overline(x)_S - c)_2^2$. Because $|S| >= 1$ and $norm(overline(x)_S - c)_2^2 >= 0$ with equality if and only if $c = overline(x)_S$, the centroid uniquely minimizes the within-cluster inertia.
]

#theorem(title: "Strict Lyapunov Descent and Finite-Step Convergence", origin: "Original Proof")[
  Let ${(cal(C)^((t)), bold(mu)^((t)))}$ be the sequence generated by Algorithm 1 with deterministic tie-breaking and non-empty clusters. Then at every iteration $t >= 0$:
  1. $J(cal(C)^((t)), bold(mu)^((t+1))) <= J(cal(C)^((t)), bold(mu)^((t))) - sum_(j=1)^k |C_j^((t))| norm(mu_j^((t+1)) - mu_j^((t)))_2^2$.
  2. $J(cal(C)^((t+1)), bold(mu)^((t+1))) <= J(cal(C)^((t)), bold(mu)^((t+1)))$, with strict inequality whenever any point switches cluster to a strictly closer centroid.
  3. Algorithm 1 terminates in at most $T <= S(n, k) <= k^n / k!$ iterations at a discrete centroidal Voronoi equilibrium, where $S(n, k)$ is the Stirling number of the second kind.
]

#proof[
  At iteration $t$, holding $cal(C)^((t))$ fixed and applying Lemma 1 independently to each cluster $C_j^((t))$ with $c = mu_j^((t))$ gives
  $ J(cal(C)^((t)), bold(mu)^((t))) - J(cal(C)^((t)), bold(mu)^((t+1))) = sum_(j=1)^k |C_j^((t))| norm(mu_j^((t+1)) - mu_j^((t)))_2^2 >= 0, $
  establishing Claim 1 with exact equality. Next, during the assignment step at iteration $t+1$, each point $x_i$ is reassigned to $c_i^((t+1)) = arg min_ell norm(x_i - mu_ell^((t+1)))_2^2$, hence
  $ norm(x_i - mu_(c_i^((t+1)))^((t+1)))_2^2 <= norm(x_i - mu_(c_i^((t)))^((t+1)))_2^2 space forall i in {1, dots, n}. $
  Summing over $i = 1, dots, n$ proves Claim 2. Combining Claims 1 and 2, the sequence $J_t := J(cal(C)^((t)), bold(mu)^((t+1)))$ is monotonically non-increasing.

  To rule out cycling, suppose $cal(C)^((t+1)) != cal(C)^((t))$ while $J_(t+1) = J_t$. Then Claim 1 implies $bold(mu)^((t+2)) = bold(mu)^((t+1))$, and Claim 2 implies every reassigned index $i$ had equal distance to its old and new centroid. Because tie-breaking deterministically selects the *smallest* index $j$ attaining the minimum distance, $c_i^((t+1)) = c_i^((t))$ once centroids are stationary, contradicting $cal(C)^((t+1)) != cal(C)^((t))$. Therefore every visited partition prior to termination corresponds to a strictly lower value of $J^*(cal(C))$ and cannot recur. Since there are only finitely many partitions of $n$ points into at most $k$ non-empty subsets—bounded by $sum_(m=1)^k S(n, m) <= k^n / k!$ for $k >= 2$—the algorithm must halt in finitely many steps at a discrete CVT.
]

#remark(title: "Bregman Generalization")[
  Banerjee et al. @banerjee2005clustering proved that Lemma 1 generalizes verbatim to *any* Bregman divergence $D_phi(x || mu)$: the arithmetic mean uniquely minimizes expected $D_phi$-distortion if and only if the distortion is a Bregman divergence. Thus Theorem 1 extends immediately to KL-divergence, Itakura–Saito, and Mahalanobis $k$-means.
]

#subsec("B", "Gauss–Newton Interpretation and Local Quadratic Rate")

#theorem(title: "Exact Block-Newton Equivalence of Lloyd Updates", origin: "Classical / Adapted — Bottou & Bengio (1995)")[
  Fix a partition $cal(C)$ in an open neighborhood $cal(U)$ of a strict CVT $bold(mu)^*$ where no data point lies on a Voronoi boundary hyperplane, i.e., $norm(x_i - mu_(c_i)^*)_2 < norm(x_i - mu_ell^*)_2$ for all $ell != c_i$. On $cal(U)$, the reduced objective $F(bold(mu)) := min_(cal(C)) J(cal(C), bold(mu))$ is $C^infinity$-smooth, and a single Lloyd update coincides identically with a full Newton–Raphson step @bottou1995convergence:
  $ bold(mu)^((t+1)) = bold(mu)^((t)) - (nabla^2 F(bold(mu)^((t))))^(-1) nabla F(bold(mu)^((t))). $
]

#proof[
  Because strict Voronoi inequalities hold on $cal(U)$ and $cal(X)$ is finite, the optimal assignment $c_i(bold(mu)) equiv c_i^*$ is locally constant on $cal(U)$. Thus $F(bold(mu)) = sum_(j=1)^k sum_(i in C_j^*) norm(x_i - mu_j)_2^2$ is a strictly convex separable quadratic form on $cal(U)$. Differentiating with respect to block $mu_j in bb(R)^d$ yields:
  $ nabla_(mu_j) F(bold(mu)) = 2 sum_(i in C_j^*) (mu_j - x_i) = 2 |C_j^*| (mu_j - overline(x)_(C_j^*)), $
  and the Hessian is strictly block-diagonal with scalar-multiple identity blocks:
  $ nabla^2_(mu_j, mu_ell) F(bold(mu)) = 2 |C_j^*| delta_(j ell) I_d. $
  Substituting into the Newton update gives
  $ mu_j^((t+1)) = mu_j^((t)) - (2 |C_j^*| I_d)^(-1) [2 |C_j^*| (mu_j^((t)) - overline(x)_(C_j^*))] = overline(x)_(C_j^*), $
  which is precisely Lloyd's update step and converges in *one single step* once $bold(mu)^((t))$ enters the interior Voronoi basin $cal(U)$. More generally, for continuous population densities $p(x)$ where Voronoi boundaries shift smoothly with $bold(mu)$, Bottou and Bengio @bottou1995convergence show this yields superlinear-to-quadratic local convergence governed by the boundary curvature integral.
]

#subsec("C", "Global NP-Hardness, Iteration Bounds & Unbounded Suboptimality")

Despite finite termination and fast local convergence, Lloyd's algorithm offers *no* global optimality guarantee.

#theorem(title: "Unbounded Approximation Ratio of Lloyd's Algorithm", origin: "Original Counterexample Proof")[
  For any target factor $M > 0$ and any $k >= 3$, there exists a one-dimensional dataset $cal(X) subset.eq bb(R)$ of size $n = 4$ and an initialization $bold(mu)^((0))$ consisting of actual data points (Forgy seeding) such that Lloyd's algorithm converges to a stationary CVT $(cal(C)_"Lloyd", bold(mu)_"Lloyd")$ satisfying
  $ J(cal(C)_"Lloyd", bold(mu)_"Lloyd") / J(cal(C)^*, bold(mu)^*) > M. $
]

#proof[
  Set $k = 3$ and construct the four-point one-dimensional dataset
  $ x_1 = 0, quad x_2 = 1, quad x_3 = L, quad x_4 = L + D, $
  with positive spacing parameters $D > 0$ and $L > 2 D + 2$. Initialize $bold(mu)^((0)) = (x_1, x_2, x_3) = (0, 1, L)$, which is a valid Forgy initialization consisting entirely of dataset elements. Then:
  - $x_1 = 0$ is nearest to $mu_1 = 0$; $x_2 = 1$ is nearest to $mu_2 = 1$; and both $x_3 = L$ and $x_4 = L + D$ are nearest to $mu_3 = L$, since $|L - 1| = L - 1 > D$.
  - After the mean update, $mu_1^((1)) = 0$, $mu_2^((1)) = 1$, and $mu_3^((1)) = L + D / 2$.
  - Check Voronoi stability for all four points under $bold(mu)^((1))$:
    - $|x_1 - mu_1^((1))| = 0 < 1$, $|x_2 - mu_2^((1))| = 0 < 1$,
    - $|x_3 - mu_3^((1))| = D / 2 < L - 1 = |x_3 - mu_2^((1))|$ (since $L > 2D + 2$),
    - $|x_4 - mu_3^((1))| = D / 2 < L + D - 1 = |x_4 - mu_2^((1))|$.
  Hence $cal(C)_"Lloyd" = {{x_1}, {x_2}, {x_3, x_4}}$ is a *strict fixed point* of Lloyd's algorithm with cost
  $ J(cal(C)_"Lloyd") = 0 + 0 + (D / 2)^2 + (D / 2)^2 = D^2 / 2. $
  In contrast, the optimal partition $cal(C)^* = {{x_1, x_2}, {x_3}, {x_4}}$ attains
  $ J(cal(C)^*) = (1 / 2)^2 + (1 / 2)^2 + 0 + 0 = 1 / 2. $
  Therefore,
  $ J(cal(C)_"Lloyd") / J(cal(C)^*) = (D^2 / 2) / (1 / 2) = D^2. $
  Choosing $D = sqrt(M) + 1$ and $L = 2D + 3$ yields $J(cal(C)_"Lloyd") / J(cal(C)^*) > M$ for arbitrarily large $M$. For general $k > 3$, append $k - 3$ isolated points at $10 L, 20 L, dots$ as singleton clusters initialized on themselves.
]

#theorem(title: "NP-Hardness, Exponential Iterations & Smoothed Runtime", origin: "Classical / Cited")[
  1. *NP-Hardness:* Minimizing $J(cal(C), bold(mu))$ is NP-hard even for $k = 2$ in general dimension $d$ @aloise2009np and even in the Euclidean plane $d = 2$ for arbitrary $k$ @mahajan2012planar.
  2. *Exponential Iteration Lower Bound:* There exist planar instances ($d = 2$) on which Lloyd's algorithm requires $2^(Omega(sqrt(n)))$ iterations before terminating @vattani2011kmeans.
  3. *Smoothed Polynomial Complexity:* If each point $x_i in [0,1]^d$ is perturbed by independent isotropic Gaussian noise $cal(N)(0, sigma^2 I_d)$, the expected number of iterations of Lloyd's algorithm is bounded by $"poly"(n, k, d, 1/sigma)$ @arthur2011smoothed.
]

#proof(label: "Proof Synthesis (Cited)")[
  Aloise et al. @aloise2009np reduce from *Densest Cut / Graph Partitioning* by embedding incidence vectors of a graph into $bb(R)^(|E|)$ so that 2-means cost corresponds directly to balanced cut maximization. Mahajan et al. @mahajan2012planar reduce from planar 3-SAT via geometric variable/clause gadgets in $bb(R)^2$. Vattani @vattani2011kmeans constructs a cascading binary-counter gadget in the plane where shifting a single centroid triggers an exponential chain of Voronoi flips. Conversely, Arthur, Manthey, and Röglin @arthur2011smoothed prove via smoothed analysis that Gaussian perturbation bounds the probability that any Voronoi bisector lies within distance $epsilon$ of a point, ensuring a minimum expected objective drop per non-trivial iteration and yielding polynomial smoothed runtime.
]

#subsec("D", "Geometric Invariances and Breakdown Robustness")

#proposition(title: "Isometry Invariance, Scale Sensitivity, and Breakdown Point 1/n", origin: "Original Proof")[
  Let $T(x) = Q x + b$ be any rigid Euclidean isometry ($Q in O(d)$, $b in bb(R)^d$), and let $S(x) = D x$ be an anisotropic diagonal scaling with $D = "diag"(s_1, dots, s_d)$, $s_r != s_q > 0$. Then:
  1. The set of stationary CVTs and global minimizers of $J$ is *strictly invariant* under $T$.
  2. Optimal partitions of $J$ are *not invariant* under $S$.
  3. The finite-sample Donoho–Huber breakdown point of the $k$-means centroids is $epsilon^* = 1 / n$.
]

#proof[
  (1) Because $Q^T Q = I_d$, for every $x_i, mu_j$ we have
  $ norm((Q x_i + b) - (Q mu_j + b))_2^2 = norm(Q (x_i - mu_j))_2^2 = norm(x_i - mu_j)_2^2. $
  Moreover, the arithmetic mean commutes with affine maps: $overline(T(x))_(C_j) = Q overline(x)_(C_j) + b$. Hence both the Voronoi step and update step commute identically with $T$.

  (2) Take four points in $bb(R)^2$: $(0, 0), (0, 2), (3, 0), (3, 2)$ with $k = 2$. Under identity scaling, vertical pairs minimize $J$ ($J = 4$ vs. horizontal pairs $J = 9$). Applying $D = "diag"(1, 10)$ stretches vertical separation to $20$, making horizontal pairs strictly optimal ($J = 9$ vs. $J = 400$). Hence anisotropic feature scaling alters the clustering.

  (3) Corrupt a single point $x_1$ by sending $norm(x_1)_2 -> infinity$ while fixing $x_2, dots, x_n$. Suppose $x_1 in C_j$ with $|C_j| = m <= n$. Then
  $ norm(mu_j)_2 = norm(1 / m x_1 + 1 / m sum_(i in C_j backslash {1}) x_i)_2 >= 1 / m norm(x_1)_2 - M_0 -> infinity, $
  where $M_0 = max_(i >= 2) norm(x_i)_2 < infinity$. Thus altering a single fraction $1/n$ of the sample drives at least one centroid outside every bounded ball in $bb(R)^d$.
]

#subsec("E", "Strong Statistical Consistency")

#theorem(title: "Strong Consistency of Empirical K-Means", origin: "Classical / Cited — Pollard (1981)")[
  Let $X_1, X_2, dots$ be i.i.d. draws from a probability distribution $P$ on $bb(R)^d$ with $bb(E)[norm(X)_2^2] < infinity$. Suppose the population objective $Phi(bold(mu)) = bb(E)[min_(1 <= j <= k) norm(X - mu_j)_2^2]$ has a unique minimizer $bold(mu)^*$ up to permutation of coordinates. Then the empirical global minimizer $hat(bold(mu))_n$ satisfies $hat(bold(mu))_n -> bold(mu)^*$ almost surely (up to permutation) as $n -> infinity$ @pollard1981strong.
]

#proof(label: "Proof Outline (Cited — Pollard 1981)")[
  Pollard @pollard1981strong first proves via the finite second moment condition that all empirical minimizers $hat(bold(mu))_n$ eventually lie inside a compact ball $B_R(0)^k$ almost surely (if any centroid diverged to infinity, relocating it to the support of $P$ strictly lowers empirical WCSS). On $B_R(0)^k$, the class of functions $g_(bold(mu))(x) = min_j norm(x - mu_j)_2^2$ is a parametric Glivenko–Cantelli class with integrable envelope $2(norm(x)_2^2 + R^2)$, implying uniform almost-sure convergence $sup_(bold(mu) in B_R^k) |hat(Phi)_n(bold(mu)) - Phi(bold(mu))| -> 0$. Standard M-estimation continuity with uniqueness of $bold(mu)^*$ yields $hat(bold(mu))_n -> bold(mu)^*$ a.s.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of Algorithm 1 evaluates $n k$ squared Euclidean distances in $bb(R)^d$ costing $Theta(n k d)$ floating-point operations and $O(n d)$ centroid accumulations. Memory consumption is $O((n + k) d)$—strictly linear in $n$, making Lloyd's algorithm vastly more memory-efficient than $O(n^2)$ distance-matrix methods (hierarchical, spectral, PAM). Because Euclidean distance squares are expanded as $norm(x_i - mu_j)_2^2 = norm(x_i)_2^2 + norm(mu_j)_2^2 - 2 x_i^T mu_j$, modern implementations formulate the assignment step as a dense BLAS Level-3 matrix multiplication $X M^T$, achieving near-peak SIMD/GPU throughput.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Lloyd's K-Means],
    [Time / Space Complexity], [$O(n k d T)$ time; $O((n+k)d)$ space; worst-case $T = 2^(Omega(sqrt(n)))$, smoothed $T = "poly"(n,k,d,1/sigma)$.],
    [Geometric Inductive Bias], [Convex, isotropic spherical Voronoi polyhedra of comparable variance.],
    [Determinism / Seeding], [Deterministic given $bold(mu)^((0))$; highly sensitive to random seeding (unbounded ratio, Theorem 3).],
    [Invariances], [Invariant under orthogonal $Q in O(d)$ & translation; sensitive to feature scaling.],
    [Outlier Breakdown Point], [$epsilon^* = 1/n$ (non-robust; single extreme outlier corrupts centroid).],
    [Recommended Upgrades], [K-Means++ (Alg. 2) for $O(log k)$ seeding; Elkan/Hamerly (Algs. 13–14) for exact speedup; PAM (Alg. 4) or Trimmed (Alg. 26) for outliers.],
  ),
  caption: [Summary of theoretical guarantees and operational properties of Lloyd's K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Because Lloyd's algorithm optimizes minimum mean-square distortion (MMSE) with $O(n k d)$ streaming-friendly linear algebra, it remains one of the highest-impact numerical primitives across engineering and science—provided the target clusters respect approximately convex, balanced, scaled Euclidean geometry.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Signal Compression & Pulse-Code Modulation (PCM) Quantization:* Originally invented by Stuart Lloyd at Bell Labs (1957, published @lloyd1982least) for optimal scalar/vector quantization. Minimizing $J$ directly minimizes expected signal-to-quantization-noise ratio (SQNR) distortion $bb(E)norm(X - Q(X))_2^2$, directly motivating the LBG vector quantizer (Algorithm 30).
2. *Billion-Scale Vector Search (IVF-PQ Indexing):* Modern retrieval-augmented generation (RAG) and dense embedding databases (e.g., FAISS, ScaNN) partition billion-vector collections into $k approx sqrt(n)$ coarse Voronoi cells via Lloyd/K-Means++, followed by subspace Product Quantization (PQ). Because Voronoi cells are convex half-space intersections, pruning search to the top-$m$ closest centroids cuts query latency by $10^2$–$10^3 times$ with high recall.
3. *Computer Vision & Color/Texture Codebooks:* Image color palette quantization (reducing 24-bit RGB to $k$ palette entries) and Bag-of-Visual-Words (BoVW) SIFT/ORB descriptor dictionaries rely on K-Means because perceptual RGB/CIELAB quantization error closely tracks squared $L_2$ distortion.
4. *Enterprise Customer & Telemetry Segmentation:* After Z-score standardization (required by Proposition 1(2)) and outlier winsorization (required by Proposition 1(3)), K-Means yields immediately interpretable centroid profiles ("average spend, tenure, frequency") for marketing cohorts, energy load-curve clustering, and cloud workload tiering.
5. *Scientific Preprocessing & Multi-Stage Pipelines:* K-Means serves as the final discretization step of Spectral Clustering (Algorithms 103–105), DeepCluster (Algorithm 191), and coreset compression (BICO Alg. 78, StreamKM++ Alg. 195), where upstream nonlinear embeddings have already mapped non-convex manifolds into linearly separable spherical clusters.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Lloyd's algorithm systematically fails on (a) *non-convex manifolds* (interlocking arcs, concentric rings), because Voronoi boundaries are strictly piecewise-linear hyperplanes; (b) *strongly anisotropic or unequal-size clusters*, where equal-radius bisectors slice large clusters into fragments (resolved by GMM Alg. 81 or Gustafson–Kessel Alg. 153); (c) *heavy-tailed outlier contamination* (resolved by PAM Alg. 4, K-Medians Alg. 15, or Trimmed K-Means Alg. 26); and (d) *purely categorical/discrete spaces* where arithmetic means are undefined (resolved by K-Modes Alg. 7 or LCA Alg. 90).

#figure(
  table(
    columns: (1.05fr, 1.15fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why K-Means Fits / Prerequisite], [Caution / Alternative],
    [ANN / Vector DBs], [Convex Voronoi pruning; BLAS-3 speed], [Use Spherical (Alg. 12) for cosine],
    [PCM / Color Quant.], [Exact MMSE distortion objective], [Use LBG (Alg. 30) splitting],
    [Customer Analytics], [Interpretable centroids after Z-score], [Use PAM (Alg. 4) if outliers exist],
    [Geospatial Hotspots], [Fast coarse zoning of dense blobs], [Use DBSCAN (Alg. 31) for arcs/noise],
    [Genomics / Microarray], [Fast baseline on PCA-reduced expression], [Use Ward (Alg. 59) or GMM (Alg. 81)],
  ),
  caption: [Application suitability and domain engineering guide for Lloyd's K-Means.],
)

#sec("V", "Structural Relation to Adjacent Paradigms")

Lloyd's K-Means serves as the canonical mathematical anchor from which dozens of methods in this monograph arise via targeted relaxations or limit regimes:

#proposition(title: "Spectral Trace Reformulation & Zero-Temperature Limit", origin: "Original Unified Derivation")[
  Let $X in bb(R)^(n times d)$ have rows $x_i^T$, and define the normalized cluster-indicator matrix $H in bb(R)^(n times k)$ by $H_(i j) = |C_j|^(-1/2)$ if $i in C_j$ and $0$ otherwise (so $H^T H = I_k$). Then minimizing $J^*(cal(C))$ is strictly equivalent to the trace maximization problem
  $ min_(cal(C)) J^*(cal(C)) = norm(X)_F^2 - max_(H in cal(H)_"ind") "Tr"(H^T X X^T H). $
  Furthermore, relaxing the discrete indicator set $cal(H)_"ind"$ to arbitrary orthonormal matrices $H^T H = I_k$ yields top-$k$ Principal Component Analysis (PCA) via the Ky Fan theorem, proving that PCA provides a tight continuous spectral lower bound on optimal $k$-means inertia.
]

#proof[
  Observe that $(H H^T)_(i ell) = |C_j|^(-1)$ if $i, ell in C_j$ and $0$ otherwise. Thus $(H H^T X)_(i,: ) = overline(x)_(C_(c_i))^T$ is the matrix whose $i$-th row equals the centroid of $x_i$'s cluster. Using idempotence $(I_n - H H^T)^2 = I_n - H H^T$,
  $ J^*(cal(C)) &= norm(X - H H^T X)_F^2 = "Tr"(X^T (I_n - H H^T) X) \
                &= norm(X)_F^2 - "Tr"(H^T X X^T H). $
  Since $cal(H)_"ind" subset.eq {H in bb(R)^(n times k) : H^T H = I_k}$, maximizing over the larger Stiefel manifold yields $sum_(r=1)^k lambda_r(X X^T)$ (the top-$k$ eigenvalues of the Gram matrix $X X^T$), which therefore lower-bounds $norm(X)_F^2 - J^*(cal(C))$ for every partition $cal(C)$.
]

This trace identity directly connects Algorithm 1 to Kernel K-Means (Algorithm 16, replacing $X X^T$ with Gram matrix $K$), Spectral Clustering (Algorithms 103–105, replacing $X X^T$ with normalized graph affinities), and NMF (Algorithm 92, relaxing orthogonality to non-negativity). Simultaneously, as proven in Part IX (Comparison Chapter II), isotropic Gaussian Mixture Models (Algorithm 81) with covariance $sigma^2 I_d$ recover Lloyd's assignment rule identically in the zero-temperature limit $sigma -> 0^+$, and Fuzzy C-Means (Algorithm 151) reduces strictly to Lloyd's algorithm as the fuzzifier $m -> 1^+$.

In empirical benchmarks on the ClusterOverview 4-dataset suite (detailed in Part IX, Chapter XIII), standard Lloyd's K-Means achieves strong speed ($< 0.4$ ms on $n <= 240$) and high Silhouette on convex Iris data, yet degrades sharply on non-convex arcs and outlier-contaminated subspaces—precisely validating Theorems 3 and Proposition 1 in practice.
