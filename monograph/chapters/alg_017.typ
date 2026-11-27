#import "../template/ieee.typ": *

#alg_entry(
  id: 17,
  name: "Global K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2003",
  paradigm: "Deterministic Incremental Greedy Seeding",
  complexity: $O(n^2 k^2 d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ be a finite dataset with squared Euclidean distortion, and let $k in {1, dots, n}$ be the target cluster cardinality. Standard Lloyd's $k$-means (Alg. 1) suffers from acute initialization sensitivity and unbounded worst-case approximation error depending on the random seed $bold(mu)^((0))$ @lloyd1982least. To eliminate stochasticity and systematically explore the combinatorial partition landscape, Likas, Vlassis, and Verbeek @likas2003global introduced *Global K-Means (GKM)*, a deterministic incremental algorithm that solves the entire nested sequence of $m$-means problems for $m = 1, 2, dots, k$.

For any stage $m in {1, dots, k}$, the $m$-means within-cluster sum of squares (WCSS) functional is
$ J_m (cal(C), bold(mu)) = sum_(j=1)^m sum_(i in C_j) norm(x_i - mu_j)_2^2, $
where $cal(C) = {C_1, dots, C_m}$ partitions ${1, dots, n}$ and $bold(mu) = (mu_1, dots, mu_m) in (bb(R)^d)^m$. Given a converged $(m - 1)$-centroid solution $bold(mu)_(m-1)^* = (mu_(m-1, 1)^*, dots, mu_(m-1, m-1)^*)$ with total WCSS $E_(m-1)$, define each point's current squared distance to its nearest existing centroid:
$ d_(m-1)^2 (x_r) := min_(1 <= j <= m-1) norm(x_r - mu_(m-1, j)^*)_2^2, $
so that $E_(m-1) = sum_(r=1)^n d_(m-1)^2 (x_r)$. When a candidate $m$-th centroid is inserted at $y in bb(R)^d$, its *one-step Voronoi distortion gain* @likas2003global is defined by
$ b_m (y) := sum_(r=1)^n max(d_(m-1)^2 (x_r) - norm(x_r - y)_2^2, space 0). $

#definition(title: "Incremental Global K-Means Trajectory")[
  Let $cal(L)_m (bold(mu)^((0)))$ denote the terminal centroid vector produced by deterministic Lloyd's $m$-means (Alg. 1) initialized at $bold(mu)^((0)) in (bb(R)^d)^m$, and let $E_m (y) := J_m (cal(C)_y, cal(L)_m (bold(mu)_(m-1)^* union {y}))$ be its resulting WCSS. The *Global K-Means trajectory* ${(cal(C)_m^*, bold(mu)_m^*)}_(m=1)^k$ is defined inductively by $bold(mu)_1^* = (overline(x))$ (the unique global $1$-means centroid $overline(x) = n^(-1) sum_(r=1)^n x_r$) and, for $m = 2, dots, k$,
  $ i_m^* &= arg min_(i in {1, dots, n}) E_m (x_i), \
    bold(mu)_m^* &= cal(L)_m (bold(mu)_(m-1)^* union {x_(i_m^*)}), $
  with $E_m := E_m (x_(i_m^*))$. In the *Fast Global K-Means (FGKM)* variant @likas2003global, exhaustive evaluation of $E_m (x_i)$ is replaced by the closed-form surrogate seed $hat(i)_m = arg max_(i in {1, dots, n}) b_m (x_i)$, invoking $cal(L)_m$ only once per stage.
]

Because GKM constructs solutions sequentially from $m = 1$ to $m = k$, it produces a complete nested library of $1$- through $k$-clusterings in a single deterministic run without random restarts (Algorithm 1).

#algbox(
  title: "Global K-Means and Fast Surrogate (GKM / FGKM)",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, target $k >= 1$, mode $in {"GKM", "FGKM"}$, max Lloyd iterations $T$.],
  outputs: [Nested solutions ${(cal(C)_m^*, bold(mu)_m^*, E_m)}_(m=1)^k$.],
)[
  1. *Base stage ($m = 1$):* Set $mu_(1,1)^* = n^(-1) sum_(r=1)^n x_r$, $cal(C)_1^* = {{1, dots, n}}$, $E_1 = sum_(r=1)^n norm(x_r - mu_(1,1)^*)_2^2$.
  2. *If* mode = `"FGKM"`, precompute $D_(r i) = norm(x_r - x_i)_2^2$ for all $1 <= r, i <= n$.
  3. *For* $m = 2, dots, k$ *do*:
  4. #h(1em) Compute $d_(m-1)^2 (x_r) = min_(1 <= j <= m-1) norm(x_r - mu_(m-1, j)^*)_2^2$ for $r = 1, dots, n$.
  5. #h(1em) *If* mode = `"GKM"`:
  6. #h(2em) *For* each candidate $i in {1, dots, n}$ *do*:
  7. #h(3em) Run deterministic Lloyd $m$-means from $bold(mu)_(m-1)^* union {x_i}$ for up to $T$ steps to get $(cal(C)^((i)), bold(mu)^((i)), E^((i)))$.
  8. #h(2em) Select $i_m^* = arg min_i E^((i))$ (smallest index on tie); set $(cal(C)_m^*, bold(mu)_m^*, E_m) = (cal(C)^((i_m^*)), bold(mu)^((i_m^*)), E^((i_m^*)))$.
  9. #h(1em) *Else* (`"FGKM"`):
  10. #h(2em) Compute $b_m (x_i) = sum_(r=1)^n max(d_(m-1)^2 (x_r) - D_(r i), 0)$ for $i = 1, dots, n$.
  11. #h(2em) Set $hat(i)_m = arg max_i b_m (x_i)$ and run Lloyd $m$-means once from $bold(mu)_(m-1)^* union {x_(hat(i)_m)}$ to obtain $(cal(C)_m^*, bold(mu)_m^*, E_m)$.
  12. *Return* ${(cal(C)_m^*, bold(mu)_m^*, E_m)}_(m=1)^k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results governing Global K-Means: (1) an exact distortion-drop identity for the one-step gain functional $b_m$; (2) strict stage-wise energy descent together with a global additive-multiplicative contraction bound relative to the true $k$-means optimum $E_k^"opt"$; (3) an explicit four-point 1D counterexample proving that despite its name, Global K-Means fails to reach global optimality with ratio approaching $4$; (4) complete determinism, rigid isometry invariance, and an exact outlier budget-capture theorem; and (5) a comparative synthesis against NP-hardness, randomized $D^2$-seeding, and local swap search.

#subsec("A", "One-Step Distortion Gain and Global Contraction Bound")

#lemma(title: "Exact One-Step Distortion Drop Identity", origin: "Original Derivation")[
  Let $bold(mu)_(m-1)^*$ be any $(m - 1)$-centroid configuration with WCSS $E_(m-1) = sum_(r=1)^n d_(m-1)^2 (x_r)$. For any candidate seed $y in bb(R)^d$, let $cal(C)_y^((0))$ denote the Voronoi partition induced by $bold(mu)_y^((0)) := bold(mu)_(m-1)^* union {y}$. Then:
  1. $J_m (cal(C)_y^((0)), bold(mu)_y^((0))) = E_(m-1) - b_m (y)$ identically.
  2. The post-Lloyd objective satisfies $E_m (y) <= E_(m-1) - b_m (y)$.
  3. For every data point $x_i in cal(X)$, $b_m (x_i) >= d_(m-1)^2 (x_i) >= 0$.
]

#proof[
  For each point $x_r in cal(X)$, its squared distance to $bold(mu)_y^((0))$ is $min(d_(m-1)^2 (x_r), norm(x_r - y)_2^2)$. Applying $min(a, c) = a - max(a - c, 0)$ pointwise gives
  $ &J_m (cal(C)_y^((0)), bold(mu)_y^((0))) = sum_(r=1)^n d_(m-1)^2 (x_r) \
    &quad - sum_(r=1)^n max(d_(m-1)^2 (x_r) - norm(x_r - y)_2^2, 0) = E_(m-1) - b_m (y), $
  establishing Claim 1. Because Lloyd's algorithm (Alg. 1) monotonically decreases $J_m$ at every assignment and centroid update step, the converged solution satisfies $E_m (y) <= J_m (cal(C)_y^((0)), bold(mu)_y^((0))) = E_(m-1) - b_m (y)$, proving Claim 2. For Claim 3, set $y = x_i$. Every summand in $b_m (x_i)$ is non-negative, while the $r = i$ summand equals $max(d_(m-1)^2 (x_i) - 0, 0) = d_(m-1)^2 (x_i)$. Dropping all $r != i$ terms yields $b_m (x_i) >= d_(m-1)^2 (x_i) >= 0$.
]

#theorem(title: "Strict Descent and Bicriteria Contraction Bound", origin: "Original Proof")[
  Let $E_k^"opt"$ be the global minimum $k$-means WCSS on $cal(X)$ (with at least $k$ distinct points). Both exact Global K-Means (GKM) and Fast Global K-Means (FGKM) satisfy, for every stage $m >= 2$:
  1. *Strict Monotonicity:* Whenever $E_(m-1) > 0$,
     $ E_m <= E_(m-1) - max_i b_m (x_i) <= (1 - 1 / n) E_(m-1) < E_(m-1). $
  2. *Global Contraction Relative to $E_k^"opt"$:*
     $ E_m <= (1 - 1 / k) E_(m-1) + 2 / k E_k^"opt". $
  3. *Closed-Form Stage-$m$ Envelope:* For any $m >= 1$ (including oversampling $m = k + s$ with $s >= 0$),
     $ E_m <= 2 E_k^"opt" + (1 - 1 / k)^(m-1) max(E_1 - 2 E_k^"opt", 0). $
]

#proof[
  At stage $m$, both GKM and FGKM satisfy $E_m <= E_m (x_(hat(i)_m)) <= E_(m-1) - max_i b_m (x_i)$ by Lemma 1(2). By Lemma 1(3) and averaging over $i = 1, dots, n$,
  $ max_i b_m (x_i) >= max_i d_(m-1)^2 (x_i) >= 1 / n sum_(r=1)^n d_(m-1)^2 (x_r) = E_(m-1) / n, $
  yielding Claim 1 immediately.

  For Claim 2, let $cal(O)^* = {O_1^*, dots, O_k^*}$ be an optimal $k$-partition of ${1, dots, n}$ with centroids $c_ell^*$ and inertias $W (O_ell^*) = sum_(r in O_ell^*) norm(x_r - c_ell^*)_2^2$, so $sum_(ell=1)^k W (O_ell^*) = E_k^"opt"$. Fix any non-empty $O_ell^*$ and any $i in O_ell^*$. Dropping non-negative terms $r in.not O_ell^*$ from $b_m (x_i)$ and using $max(a, 0) >= a$ gives $b_m (x_i) >= sum_(r in O_ell^*) (d_(m-1)^2 (x_r) - norm(x_r - x_i)_2^2)$. Averaging uniformly over $i in O_ell^*$ yields
  $ max_i b_m (x_i) &>= sum_(r in O_ell^*) d_(m-1)^2 (x_r) \
    &quad - 1 / (|O_ell^*|) sum_(i, r in O_ell^*) norm(x_r - x_i)_2^2. $
  By the pairwise inertia identity (inserting $plus.minus c_ell^*$ with $sum_(r in O_ell^*) (x_r - c_ell^*) = bold(0)$), $sum_(i, r in O_ell^*) norm(x_r - x_i)_2^2 = 2 |O_ell^*| W (O_ell^*)$, hence $max_i b_m (x_i) >= sum_(r in O_ell^*) d_(m-1)^2 (x_r) - 2 W (O_ell^*)$. Summing over $ell = 1, dots, k$ gives
  $ k max_i b_m (x_i) >= sum_(r=1)^n d_(m-1)^2 (x_r) - 2 sum_(ell=1)^k W (O_ell^*) = E_(m-1) - 2 E_k^"opt". $
  Dividing by $k$ and combining with $E_m <= E_(m-1) - max_i b_m (x_i)$ establishes Claim 2:
  $ E_m <= E_(m-1) - (E_(m-1) - 2 E_k^"opt") / k = (1 - 1 / k) E_(m-1) + 2 / k E_k^"opt". $

  Subtracting $2 E_k^"opt"$ from both sides produces the linear contraction
  $ E_m - 2 E_k^"opt" <= (1 - 1 / k) (E_(m-1) - 2 E_k^"opt"). $
  Because $1 - 1 / k >= 0$, induction from $m = 1$ gives $E_m - 2 E_k^"opt" <= (1 - 1 / k)^(m-1) (E_1 - 2 E_k^"opt") <= (1 - 1 / k)^(m-1) max(E_1 - 2 E_k^"opt", 0)$, completing the proof of Claim 3.
]

#remark(title: "Bicriteria Oversampling Consequence")[
  Because $(1 - 1/k)^(m-1) <= exp(-(m - 1)/k)$, Theorem 1(3) proves that running deterministic Global K-Means (or FGKM) for $m = k + ceil(k ln(E_1 / (epsilon E_k^"opt")))$ clusters achieves a guaranteed $(2 + epsilon)$-approximation to the optimal $k$-means cost—an explicit deterministic bicriteria guarantee.
]

#subsec("B", "Failure of Global Optimality: Exact 4-Point Counterexample")

Although its name suggests global optimality, incremental seeding locks the algorithm into the optimal $(m - 1)$-centroid topology before stage $m$ begins. We prove constructively that even on four points in $bb(R)$ with $k = 3$, Global K-Means fails to find the global optimum under *every* tie-breaking rule.

#theorem(title: "Suboptimality of Global K-Means on n = 4, d = 1, k = 3", origin: "Original Counterexample")[
  For any parameter $r in (3, 5)$, let $cal(X)_r = {-r, -1, +1, +r} subset.eq bb(R)$ and $k = 3$. Then, independent of tie-breaking, exact Global K-Means converges to a strictly suboptimal partition with
  $ E_3^"GKM" / E_3^"opt" = (r - 1)^2 / 4, $
  which approaches $4$ as $r -> 5^-$ (and exceeds $2.25$ at $r = 4$).
]

#proof[
  First determine $E_3^"opt"$. Any $3$-partition of four points pairs exactly two points into a cluster of size $2$ and leaves two singletons. For $r in (3, 5)$, pairing ${-1, +1}$ costs $1^2 + 1^2 = 2$, whereas pairing ${-r, -1}$ or ${+1, +r}$ costs $(r - 1)^2 / 2 > 2$ (since $r > 3$), and all non-adjacent pairs cost strictly more. Thus $cal(C)_3^"opt" = {{-r}, {-1, +1}, {+r}}$ is unique with $E_3^"opt" = 2$.

  Now trace Algorithm 1 stage by stage on $cal(X)_r$:
  - *Stage $m = 1$:* The global centroid is $mu_(1,1)^* = 0$ with $E_1 = 2 r^2 + 2$.
  - *Stage $m = 2$:* Evaluate all four candidates $y in {-r, -1, +1, +r}$ alongside $mu_(1,1)^* = 0$:
    - For $y = -r$: since $|-1 - 0| = 1 < r - 1 = |-1 - (-r)|$, Voronoi assignment yields ${-r}$ and ${-1, +1, +r}$. Updating means gives $(-r, r/3)$. Since $|-1 - r/3| = 1 + r/3 < r - 1$ for all $r > 3$, Lloyd halts immediately with WCSS $2 + 2 r^2 / 3$.
    - For $y = -1$: Voronoi assignment maps ${-r, -1} -> -1$ and ${+1, +r} -> 0$ (since $|1 - 0| = 1 < 2$). Mean update gives $(-(r+1)/2, +(r+1)/2)$, which is a strict fixed point with WCSS $4 ((r - 1)/2)^2 = (r - 1)^2$.
    - By reflection symmetry, $y = +1$ gives $(r - 1)^2$ and $y = +r$ gives $2 + 2 r^2 / 3$.
    Compare the two outcomes: $(r - 1)^2 < 2 + 2 r^2 / 3$ holds if and only if $r^2 - 6 r - 3 < 0$, i.e., $r < 3 + 2 sqrt(3) approx 6.464$. Since $r in (3, 5)$, Global K-Means uniquely selects $bold(mu)_2^* = (-(r+1)/2, +(r+1)/2)$ with $E_2 = (r - 1)^2$.
  - *Stage $m = 3$:* Starting from $bold(mu)_2^*$, test all four seeds:
    - For $y = -r$: initial centers are $(-(r+1)/2, +(r+1)/2, -r)$. Since $|-1 - (-(r+1)/2)| = (r - 1)/2 < r - 1$, Voronoi assignment gives ${-r}$, ${-1}$, ${+1, +r}$. Mean update produces $(-1, +(r+1)/2, -r)$. Now check $+1$: its distance to $+(r+1)/2$ is $(r - 1)/2$, while its distance to $-1$ is $2$. Because $r < 5$, $(r - 1)/2 < 2$ *strictly*, so $+1$ remains in ${+1, +r}$ and Lloyd halts with WCSS $(r - 1)^2 / 2$.
    - For $y = -1$: Voronoi assignment assigns $-r -> -(r+1)/2$, $-1 -> -1$, and $+1, +r -> +(r+1)/2$ (again because $(r - 1)/2 < 2$). Mean update gives the identical centroid set $(-r, +(r+1)/2, -1)$, halting at WCSS $(r - 1)^2 / 2$.
    - By reflection symmetry, $y = +1$ and $y = +r$ yield the reflected partition ${-r, -1}, {+1}, {+r}$ with identical WCSS $(r - 1)^2 / 2$.

  No ties occur at any assignment step. Therefore $E_3^"GKM" = (r - 1)^2 / 2$, giving $E_3^"GKM" / E_3^"opt" = (r - 1)^2 / 4 -> 4^-$ as $r -> 5^-$.
]

#subsec("C", "Invariance, Determinism, and Outlier Budget Capture")

#proposition(title: "Isometry Invariance and Outlier Budget Collapse", origin: "Original Proof")[
  Assume $cal(X) subset.eq bb(R)^d$ is in general position (no equidistant ties). Then:
  1. Both GKM and FGKM are *strictly deterministic*, invariant under arbitrary data permutations $pi in S_n$, and invariant under all rigid Euclidean isometries $T(x) = Q x + b$ ($Q in O(d)$, $b in bb(R)^d$), yet *not* invariant under anisotropic feature scaling.
  2. The Donoho–Huber centroid breakdown point is $epsilon^* = 1 / n$. Specifically, if $n - 1$ clean points lie in $B_R (0)$ and a single outlier $x_n$ satisfies $norm(x_n)_2 -> infinity$, then for $k >= 2$, GKM deterministically allocates one full centroid exclusively to the singleton ${x_n}$ at stage $m = 2$, reducing the effective cluster budget for the inliers from $k$ to $k - 1$.
]

#proof[
  (1) Because $Q^T Q = I_d$, $norm(T(x_r) - T(y))_2^2 = norm(x_r - y)_2^2$ for all $x_r, y$, and arithmetic means commute with affine maps: $overline(T(x))_C = Q overline(x)_C + b$. Consequently, every pairwise squared distance $D_(r i)$, every gain $b_m (x_i)$, and every Lloyd trajectory $cal(L)_m$ commutes identically with $T$. Because Stage $m$ minimizes over the *unordered* set $cal(X)$ and general position excludes ties, the selected seed and partition are strictly independent of sample ordering $pi$. Anisotropic scaling fails invariance via the same four-point rectangle counterexample as Alg. 1.

  (2) Let $overline(x)_"clean"$ and $E_"clean" <= (n - 1) R^2$ denote the mean and $1$-means inertia of ${x_1, dots, x_(n-1)} subset.eq B_R (0)$, and set $M := norm(x_n - overline(x)_"clean")_2 -> infinity$. At $m = 1$, $mu_(1,1)^* = overline(x)_"clean" + n^(-1) (x_n - overline(x)_"clean")$, so $norm(mu_(1,1)^*)_2 = Theta(M) -> infinity$, establishing breakdown point $1 / n$. At stage $m = 2$, testing candidate $y = x_n$ yields initial Voronoi partition ${x_n}$ and ${x_1, dots, x_(n-1)}$ once $M > 2 n R / (n - 2)$ (since for every $r < n$, $norm(x_r - mu_(1,1)^*)_2 <= 2 R + M / n < M - 2 R <= norm(x_r - x_n)_2$). One mean update produces centroids $(overline(x)_"clean", x_n)$ with total WCSS $E_"clean" = O(R^2)$, whereas any $2$-partition keeping $x_n$ grouped with at least one clean point $x_r$ incurs inertia at least $norm(x_n - x_r)_2^2 / 2 >= (M - R)^2 / 2 -> infinity$. Thus for all sufficiently large $M$, GKM strictly isolates ${x_n}$ at $m = 2$, permanently dedicating one centroid to the outlier across stages $m = 2, dots, k$.
]

#subsec("D", "Classical Context: Hardness, Local Search, and Seeding")

#theorem(title: "Complexity Landscape vs. Randomized and Swap Methods", origin: "Classical — Cited")[
  1. *NP-Hardness:* Exact minimization of $J_k$ is NP-hard for $k = 2$ in arbitrary dimension $d$ @aloise2009np and for $d = 2$ with arbitrary $k$ @mahajan2012planar.
  2. *Randomized $D^2$-Seeding:* K-Means++ (Alg. 2) samples each new center with probability proportional to $d_(m-1)^2 (x_i)$ in $O(n k d)$ time and achieves expected cost $bb(E)[E_k] <= 8 (ln k + 2) E_k^"opt"$ @arthur2007kmeanspp.
  3. *Multi-Centroid Local Search:* Swapping up to $p$ centroids simultaneously achieves a $(9 + 2 / p)$-approximation in $O(n^(p+1) k^p)$ time @kanungo2004local.
  4. *Deterministic Incremental Search:* GKM @likas2003global requires $O(n^2 k^2 d T)$ time (reduced to $O(n^2 d + n k^2 d T)$ by FGKM), eliminates RNG variance, and furnishes an entire hierarchy $m = 1, dots, k$ simultaneously.
]

#proof(label: "Proof Synthesis (Cited)")[
  Aloise et al. @aloise2009np and Mahajan et al. @mahajan2012planar rule out polynomial-time exact solvers unless $"P" = "NP"$. Arthur and Vassilvitskii @arthur2007kmeanspp prove via a potential-function induction on "covered" versus "uncovered" optimal clusters that sampling proportionally to $d_(m-1)^2$ yields an $O(log k)$ multiplicative approximation in expectation, whereas deterministic greedy insertion (Theorem 1) yields an additive-multiplicative contraction toward $2 E_k^"opt"$ without logarithmic growth in $k$ under bicriteria oversampling, yet remains susceptible to fixed-budget traps (Theorem 2). Kanungo et al. @kanungo2004local show that escaping traps like Theorem 2 with a constant multiplicative ratio requires *simultaneous* removal-and-insertion swaps rather than purely additive insertion.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Exact Global K-Means executes $n$ full runs of Lloyd's $m$-means at each stage $m = 2, dots, k$. Taking $T$ iterations per Lloyd call yields $sum_(m=2)^k n dot O(n m d T) = O(n^2 k^2 d T)$ arithmetic operations with $O((n + k) d)$ memory. While quadratic in both $n$ and $k$, the *Fast Global K-Means* surrogate reduces total work to $O(n^2 d + n k^2 d T)$ by precomputing the $n times n$ squared distance matrix once in $O(n^2 d)$ time ($O(n^2)$ memory), evaluating $b_m (x_i)$ across all $i$ in $O(n^2)$ per stage, and running Lloyd's algorithm only once per $m$.

Further large-scale speedups preserve GKM's deterministic character: (i) restricting candidates $x_i$ to $M << n$ *kd-tree* bucket representatives @likas2003global or an $epsilon$-coreset @harpeled2004coresets drops the outer search factor from $n$ to $M$; (ii) replacing inner Lloyd calls with Elkan (Alg. 13) @elkan2003using or Hamerly (Alg. 14) @hamerly2010making triangle-inequality pruning eliminates redundant distance evaluations; and (iii) because the $n$ candidate runs at stage $m$ are embarrassingly parallel, multi-core or GPU evaluation scales linearly with processor count.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Global K-Means],
    [Time / Space Complexity], [GKM: $O(n^2 k^2 d T)$ time, $O((n+k)d)$ space; FGKM: $O(n^2 d + n k^2 d T)$ time, $O(n^2 + n d)$ space.],
    [Geometric Inductive Bias], [Convex spherical Voronoi cells built via nested greedy $1 -> k$ refinement.],
    [Determinism / Seeding], [100% deterministic and permutation-invariant (no RNG seed or multi-start tuning).],
    [Approximation Profile], [Strict descent $E_m < E_(m-1)$; bicriteria $(2+epsilon)$-approx; worst-case ratio $>= 4 - epsilon$ at $k = 3$ (Thm. 2).],
    [Invariances / Robustness], [Rigid isometry invariant; scale-sensitive; $epsilon^* = 1/n$ with deterministic outlier-budget capture (Prop. 1).],
    [Recommended Upgrades], [FGKM or coreset subsampling for $n > 10^4$; K-Means++ (Alg. 2) for $O(n k d)$ streaming; PAM (Alg. 4) for outliers.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Global K-Means (GKM / FGKM).],
)

#sec("IV", "Engineering and Scientific Applications")

Global K-Means excels in regimes where sample sizes are moderate ($n <= 10^4$, or larger via FGKM/coresets), reproducibility is mandatory, and practitioners must simultaneously inspect clustering quality across every candidate cardinality $m = 1, dots, k$.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Deterministic Model-Order Selection & Elbow Curves:* Selecting $k$ via the Elbow method, Silhouette index @rousseeuw1987silhouettes, or Caliński–Harabasz score normally requires dozens of random K-Means restarts at every $m in {2, dots, k_"max"}$. Because GKM produces a strictly monotonically decreasing curve $E_1 > E_2 > dots > E_(k_"max")$ in one pass (Theorem 1), curvature diagnostics are completely free of seed-induced jitter.
2. *Regulatory & Clinical Medical Imaging Segmentation:* Radiomics, MRI tissue quantification, and histopathology pipelines audited by regulatory bodies require bit-for-bit reproducible partitions invariant to input ordering. Proposition 1(1) guarantees deterministic, rotation-invariant segmentation without random-seed artifacts.
3. *Vector Quantization Codebook Initialization:* In rate-distortion coding of speech and image blocks, evaluating nested codebook sizes $m = 1, 2, dots, k$ directly traces the empirical operational rate-distortion frontier $D(R) = E_m / n$, complementing LBG binary splitting (Alg. 30) @linde1980algorithm when non-power-of-two codebook sizes are required.
4. *Materials Science & Experimental Design:* High-throughput combinatorial chemistry and alloy screening typically analyze $n in [10^2, 5 times 10^3]$ standardized descriptors where evaluation cost is dwarfed by laboratory cost; exhaustive deterministic search over every sample as candidate prototype reliably uncovers small, compact phase clusters that random seeding frequently misses.
5. *Benchmark Oracle & Reproducible Baseline:* In clustering research and automated ML benchmarking, GKM serves as a zero-hyperparameter deterministic upper-quality reference for Euclidean WCSS against which faster stochastic algorithms (K-Means++, Mini-Batch, K-Means$||$ @bahmani2012scalable) are calibrated.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Global K-Means fails or becomes impractical under four conditions: (a) *massive datasets* ($n >= 10^5$), where exact $O(n^2 k^2 d T)$ scaling is prohibitive without coreset reduction or K-Means++ (Alg. 2); (b) *non-convex manifolds* (interlocking moons, concentric rings), inherited from Euclidean Voronoi geometry; (c) *heavy outlier contamination*, where Proposition 1(2) shows GKM greedily sacrifices cluster slots to isolate extreme singletons instead of partitioning the inlier bulk; and (d) *high-dimensional noise subspaces*, where pairwise distance concentration @beyer1999nearest flattens $b_m (x_i)$ contrasts across candidates.

#figure(
  table(
    columns: (1.05fr, 1.15fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Global K-Means Fits], [Caution / Alternative],
    [Model Selection ($k$)], [Free monotonic $E_1 > dots > E_k$ sweep], [Use FGKM if $n > 2000$],
    [Clinical / Audit ML], [Strict determinism & order invariance], [Standardize; winsorize outliers],
    [Materials / Chem.], [Exhaustive search finds small phases], [Use GMM (Alg. 81) if elliptical],
    [Web-Scale Vectors], [Coreset-GKM for coarse quantization], [Use K-Means++ (Alg. 2) natively],
    [Noisy Sensor Data], [Isolates extreme anomalies at $m = 2$], [Use Trimmed (Alg. 26) or DBSCAN],
  ),
  caption: [Application suitability and engineering trade-off matrix for Global K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Global K-Means sits at the exact mathematical intersection of stochastic $D^2$-seeding (K-Means++, Alg. 2) and deterministic top-down divisive clustering (Bisecting K-Means @steinbach2000comparison and LBG @linde1980algorithm).

#proposition(title: "Zero-Temperature Limit of D^2-Seeding and Gain Dominance", origin: "Original Unified Proof")[
  Define the inverse-temperature $D^2$-seeding distribution $p_m^((beta)) (x_i) := d_(m-1)^(2 beta) (x_i) / sum_(r=1)^n d_(m-1)^(2 beta) (x_r)$ for $beta > 0$ (recovering K-Means++ @arthur2007kmeanspp at $beta = 1$). Then:
  1. As $beta -> infinity$, $p_m^((beta))$ concentrates on $arg max_i d_(m-1)^2 (x_i)$ (the farthest-first Gonzalez seed), which maximizes only the *single-point lower bound* $d_(m-1)^2 (x_i) <= b_m (x_i)$ from Lemma 1(3).
  2. In contrast, FGKM maximizes the *full Voronoi capture integral* $b_m (x_i) = d_(m-1)^2 (x_i) + sum_(r != i) max(d_(m-1)^2 (x_r) - norm(x_r - x_i)_2^2, 0)$, simultaneously penalizing isolated outliers that capture zero neighbors while rewarding dense unrepresented modes.
]

#proof[
  Let $d_"max"^2 := max_r d_(m-1)^2 (x_r) > 0$ and let $cal(M) = {i : d_(m-1)^2 (x_i) = d_"max"^2}$. For any $i in.not cal(M)$, $(d_(m-1)^2 (x_i) / d_"max"^2)^beta -> 0$ exponentially as $beta -> infinity$, so $p_m^((beta)) (x_i) -> 0$ for $i in.not cal(M)$ and $p_m^((beta)) (x_i) -> 1 / |cal(M)|$ for $i in cal(M)$, proving Claim 1. Claim 2 follows directly by splitting the $r = i$ term out of $b_m (x_i)$: an isolated outlier $x_i$ with $norm(x_r - x_i)_2^2 >= d_(m-1)^2 (x_r)$ for all $r != i$ attains exactly $b_m (x_i) = d_(m-1)^2 (x_i)$, whereas a point inside an unrepresented cluster of $s$ tight points at distance $Delta$ from existing centroids attains $b_m (x_i) approx s Delta^2$, dominating the outlier whenever $s Delta^2 > d_(m-1)^2 (x_"out")$.
]

This structural link also clarifies GKM's relation to Hartigan–Wong (Alg. 19): GKM acts as a deterministic *global seeding outer loop* whose inner Lloyd calls can be replaced directly by Hartigan–Wong refinement to enforce stricter 1-opt equilibria at every stage $m$.

In empirical evaluation on the ClusterOverview 4-dataset benchmark suite (using Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), Global K-Means attains strong, reproducible performance across all four regimes: *Iris* $"Sil" = 0.447$, $"ARI" = 0.583$; *Customers* $"Sil" = 0.365$; *Multi-Shape* $"Sil" = 0.712$ ($"ARI" = "N/A"$); and *High-D Subspace+Noise* $"Sil" = 0.303$ ($"ARI" = "N/A"$), with a mean runtime of $0.398$ ms. Notably, while matching Hartigan–Wong (Alg. 19) identically on Iris ($0.447 / 0.583$) and Multi-Shape ($0.712$), GKM achieves meaningfully higher cluster separation on both Customers ($0.365$ vs. $0.332$) and High-D Subspace+Noise ($0.303$ vs. $0.259$ for Alg. 19 and $0.181$ for ISODATA, Alg. 18)—empirically confirming that exhaustive deterministic stage-wise search avoids poor random Voronoi basins, at a modest $8 times$ runtime premium over Hartigan–Wong ($0.398$ ms vs. $0.050$ ms) consistent with its $O(n^2)$ candidate sweep.
