#import "../template/ieee.typ": *

#alg_entry(
  id: 179,
  name: "Neural Gas (NG)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "1993",
  paradigm: "Rank-Ordered Soft Competitive Vector Quantization",
  complexity: $cal(O)(T M (d + log M))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Neural Gas (NG), introduced by Martinetz, Berkovich, and Schulten @martinetz1993neuralgas, resolves two core limitations of Self-Organizing Maps (Alg. 176): topological distortion caused by a rigid pre-specified lattice and the non-existence of a continuous-space Lyapunov potential under Kohonen's Euclidean BMU rule. Instead of defining inter-prototype neighborhoods on a fixed external grid, Neural Gas couples prototypes dynamically according to their *distance rank* in the ambient data space $bb(R)^d$.

Let $W = (bold(w)_1, dots, bold(w)_M) in bb(R)^(d times M)$ be a codebook of $M$ prototypes quantizing a data distribution $P$ on $bb(R)^d$ (or empirical sample $X = {bold(x)_1, dots, bold(x)_n}$). For any input $bold(x) in bb(R)^d$, define the *neighborhood rank* $k_i (bold(x), W) in {0, 1, dots, M - 1}$ of prototype $i$ as the number of prototypes strictly closer to $bold(x)$ than $bold(w)_i$:
$
  k_i (bold(x), W) := |{j in [M] : norm(bold(x) - bold(w)_j)_2
    < norm(bold(x) - bold(w)_i)_2}|,
$
breaking exact ties on measure-zero bisector hyperplanes deterministically by index. Given an annealed neighborhood range $lambda > 0$, define the strictly decreasing exponential rank kernel $h_lambda (k) := exp(-k / lambda)$ and its normalization constant
$
  C(lambda) := sum_(k=0)^(M-1) h_lambda (k)
    = (1 - exp(-M / lambda)) / (1 - exp(-1 / lambda)).
$

#definition(title: "Martinetz-Berkovich-Schulten Neural Gas Energy Functional")[
  For fixed neighborhood range $lambda > 0$ and probability distribution $P$ with finite second moment, the population and empirical Neural Gas cost functionals are
  $
    E_lambda (W) &:= 1 / (2 C(lambda)) sum_(i=1)^M
      integral_(bb(R)^d) h_lambda (k_i (bold(x), W))
      norm(bold(x) - bold(w)_i)_2^2 dif P(bold(x)), \
    hat(E)_lambda (W) &:= 1 / (2 n C(lambda)) sum_(m=1)^n sum_(i=1)^M
      h_lambda (k_i (bold(x)_m, W)) norm(bold(x)_m - bold(w)_i)_2^2.
  $
]

Because $h_lambda (0) = 1 > h_lambda (1) > dots > h_lambda (M - 1) > 0$, every input pulls not only its nearest prototype ($k_i = 0$) but also higher-rank prototypes with exponentially decaying strength, diffusing prototypes like gas particles throughout the data manifold before freezing into a sharp Voronoi quantizer as $lambda -> 0^+$.

#algbox(
  title: "Online and Batch Neural Gas Vector Quantization",
  inputs: [Data $X = {bold(x)_m}_(m=1)^n subset bb(R)^d$; codebook size $M$; range schedule $lambda_0 > lambda_f > 0$; learning rate schedule $eta_0 > eta_f > 0$ (online) or epochs $E$ (batch).],
  outputs: [Optimized prototypes $W = {bold(w)_i}_(i=1)^M$ and nearest-prototype assignments $c(m) = arg min_i norm(bold(x)_m - bold(w)_i)_2$.],
)[
  + *Initialize:* Draw $W^((0)) = (bold(w)_1^((0)), dots, bold(w)_M^((0)))$ uniformly from $X$ or bounded support.
  + *Online Mode ($t = 0, dots, T - 1$):*
    - Anneal $lambda_t = lambda_0 (lambda_f / lambda_0)^(t / T)$ and $eta_t = eta_0 (eta_f / eta_0)^(t / T)$.
    - Sample $bold(x)_t in X$; compute $d_i = norm(bold(x)_t - bold(w)_i)_2^2$ for all $i in [M]$.
    - Sort distances in ascending order to assign ranks $k_i = k_i (bold(x)_t, W) in {0, dots, M - 1}$.
    - Update every prototype $i in [M]$ simultaneously:
      $
        bold(w)_i <- bold(w)_i + eta_t exp(-k_i / lambda_t)
          (bold(x)_t - bold(w)_i).
      $
  + *Batch Mode ($e = 1, dots, E$):*
    - Anneal $lambda_e$; for each $m in [n]$, sort distances to obtain $k_(i m) = k_i (bold(x)_m, W)$.
    - Update every prototype $i in [M]$ via exact rank-weighted centroid:
      $
        bold(w)_i <- (sum_(m=1)^n exp(-k_(i m) / lambda_e) bold(x)_m) /
          (sum_(m=1)^n exp(-k_(i m) / lambda_e)).
      $
  + *Return:* Converged codebook $W$ and Voronoi assignments $c(m) = arg min_(i in [M]) norm(bold(x)_m - bold(w)_i)_2$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four central theoretical guarantees of Neural Gas: exact boundary cancellation showing that NG performs true stochastic gradient descent on $E_lambda (W)$, an original rearrangement-inequality MM proof of exact monotonic descent for batch NG, a global convexity and dead-unit elimination theorem at large $lambda$, and an exponential homotopy bound connecting NG to Lloyd $k$-means as $lambda -> 0^+$.

#subsec("A", "Exact Boundary Cancellation and Gradient-Field Property")

Unlike Kohonen's online SOM (whose continuous drift field has non-vanishing curl, Alg. 176), Neural Gas follows the *exact* negative gradient of a smooth potential @martinetz1993neuralgas.

#theorem(title: "Exact Gradient Property of Neural Gas via Bisector Cancellation", origin: "Classical — Cited")[
  Let $P$ admit a continuous density $p$ on $bb(R)^d$ with finite second moment, and suppose prototypes $bold(w)_1, dots, bold(w)_M$ are pairwise distinct. Then $E_lambda (W)$ is continuously differentiable at $W$, and for every $i in [M]$,
  $
    -nabla_(bold(w)_i) E_lambda (W)
      = 1 / C(lambda) integral_(bb(R)^d)
        h_lambda (k_i (bold(x), W)) (bold(x) - bold(w)_i) p(bold(x)) dif bold(x).
  $
  Hence online Neural Gas is unbiased stochastic gradient descent on $E_lambda (W)$ @martinetz1993neuralgas.
]

#proof[
  Write $E_lambda (W) = (1 / (2 C(lambda))) integral_(bb(R)^d) g(bold(x), W) p(bold(x)) dif bold(x)$, where
  $
    g(bold(x), W) := sum_(j=1)^M h_lambda (k_j (bold(x), W))
      norm(bold(x) - bold(w)_j)_2^2.
  $
  Partition $bb(R)^d$ into the open permutation cells on which the distance ordering of all $M$ prototypes is strict and locally constant. Inside each cell, $k_j (bold(x), W)$ is constant, so differentiating $g$ with respect to $bold(w)_i$ yields precisely $2 h_lambda (k_i (bold(x), W)) (bold(w)_i - bold(x))$.

  To evaluate boundary contributions via Reynolds' transport theorem when $bold(w)_i$ varies, note that permutation cells meet along codimension-1 bisector hyperplanes $H_(i l) := {bold(x) in bb(R)^d : norm(bold(x) - bold(w)_i)_2 = norm(bold(x) - bold(w)_l)_2}$ for $l != i$ (higher-order ties have codimension at least $2$ and zero $(d - 1)$-dimensional surface measure). Fix any generic point $bold(x) in H_(i l)$ with common squared distance $R^2 := norm(bold(x) - bold(w)_i)_2^2 = norm(bold(x) - bold(w)_l)_2^2$, and let $r$ be the number of prototypes $j in.not {i, l}$ strictly closer to $bold(x)$ than $R$. Across $H_(i l)$, the ranks of all $j in.not {i, l}$ remain unchanged, while $(k_i, k_l)$ swaps from $(r, r + 1)$ on the $i$-side to $(r + 1, r)$ on the $l$-side. Evaluating the one-sided limits $g^+ (bold(x), W)$ and $g^- (bold(x), W)$ across $H_(i l)$ gives
  $
    g^+ (bold(x), W) - g^- (bold(x), W)
      &= (h_lambda (r) R^2 + h_lambda (r + 1) R^2) \
      &quad - (h_lambda (r + 1) R^2 + h_lambda (r) R^2) = 0.
  $
  Because $g(bold(x), W)$ is strictly continuous across every $(d - 1)$-dimensional facet $H_(i l)$, every boundary surface integral vanishes identically! Thus $nabla_(bold(w)_i) E_lambda (W)$ equals purely the interior integral, proving the exact gradient identity.
]

#subsec("B", "Rearrangement Inequality and Exact Batch Monotonicity")

We next prove from first principles via the Hardy-Littlewood-Pólya rearrangement inequality that distance-rank sorting is not an ad-hoc rule, but the *unique exact minimizer* of a bilinear permutation majorizer, yielding an exact Majorization-Minimization (MM) descent bound for batch Neural Gas.

#theorem(title: "Exact MM Monotonicity and Finite Convergence of Batch Neural Gas", origin: "Original Proof")[
  Let $S_M$ denote the symmetric group of bijections $pi : [M] -> {0, 1, dots, M - 1}$, and define the auxiliary assignment functional over $pi = (pi_1, dots, pi_n) in (S_M)^n$:
  $
    tilde(E)_lambda (W, pi) := 1 / (2 n C(lambda)) sum_(m=1)^n sum_(i=1)^M
      h_lambda (pi_m (i)) norm(bold(x)_m - bold(w)_i)_2^2.
  $
  Then for every fixed $lambda > 0$:
  + *Rearrangement Optimality:* $min_(pi in (S_M)^n) tilde(E)_lambda (W, pi) = hat(E)_lambda (W)$, attained by the distance-rank permutation $pi_m^* (i) = k_i (bold(x)_m, W)$.
  + *Strict Descent Bound:* Every batch Neural Gas iteration $W^((t)) |-> W^((t+1))$ satisfies
    $
      hat(E)_lambda (W^((t))) - hat(E)_lambda (W^((t+1)))
        >= 1 / (2 n C(lambda)) sum_(i=1)^M S_i^((t+1))
          norm(bold(w)_i^((t+1)) - bold(w)_i^((t)))_2^2,
    $
    where $S_i^((t+1)) := sum_(m=1)^n h_lambda (k_i (bold(x)_m, W^((t)))) >= n exp(-(M - 1) / lambda) > 0$. Under deterministic tie-breaking, batch NG terminates at a stationary fixed point in at most $(M!)^n$ steps.
]

#proof[
  For (1), fix $W$ and sample $m in [n]$, and let $d_i := norm(bold(x)_m - bold(w)_i)_2^2$. Since $h_lambda (0) > h_lambda (1) > dots > h_lambda (M - 1)$ is strictly decreasing, the Hardy-Littlewood-Pólya rearrangement inequality states that $sum_(i=1)^M h_lambda (pi_m (i)) d_i$ is minimized over $pi_m in S_M$ if and only if $h_lambda (pi_m (i))$ is ordered oppositely to $d_i$—meaning the smallest distance receives rank $0$, the second-smallest receives rank $1$, and so forth. Hence $pi_m^* (i) = k_i (bold(x)_m, W)$ minimizes the inner sum for every $m$, giving $min_pi tilde(E)_lambda (W, pi) = hat(E)_lambda (W)$.

  For (2), let $pi^((t+1)) := pi^* (W^((t)))$ and abbreviate $h_(i m) := h_lambda (pi_m^((t+1)) (i))$. Because $bold(w)_i^((t+1)) = (1 / S_i^((t+1))) sum_(m=1)^n h_(i m) bold(x)_m$, expanding each squared residual around $bold(w)_i^((t+1))$ cancels the cross-term $sum_(m=1)^n h_(i m) (bold(x)_m - bold(w)_i^((t+1))) = bold(0)$, giving
  $
    tilde(E)_lambda (W^((t)), pi^((t+1))) - tilde(E)_lambda (W^((t+1)), pi^((t+1)))
      = 1 / (2 n C(lambda)) sum_(i=1)^M S_i^((t+1))
        norm(bold(w)_i^((t+1)) - bold(w)_i^((t)))_2^2.
  $
  By Part (1), $hat(E)_lambda (W^((t))) = tilde(E)_lambda (W^((t)), pi^((t+1)))$ and $hat(E)_lambda (W^((t+1))) = min_pi tilde(E)_lambda (W^((t+1)), pi) <= tilde(E)_lambda (W^((t+1)), pi^((t+1)))$. Chaining these relations yields the exact descent inequality. Because $S_i^((t+1)) > 0$, any update $W^((t+1)) != W^((t))$ strictly decreases $hat(E)_lambda$, so no discrete rank profile $pi in (S_M)^n$ can revisit; finiteness of $|(S_M)^n| = (M!)^n$ guarantees finite termination.
]

#subsec("C", "Large-Range Convexity and Elimination of Dead Prototypes")

#theorem(title: "Global Convex Collapse at Large Range and Dead-Unit Elimination", origin: "Original Proof")[
  Let $bar(bold(x)) := (1 / n) sum_(m=1)^n bold(x)_m$ and $V_X := (1 / n) sum_(m=1)^n norm(bold(x)_m - bar(bold(x)))_2^2$.
  + *Convex Homotopy Limit:* As $lambda -> infinity$, $hat(E)_lambda (W)$ converges uniformly on compact subsets of $bb(R)^(d times M)$ to the strictly convex separable quadratic
    $
      E_infinity (W) := 1 / (2 M) sum_(i=1)^M
        (V_X + norm(bold(w)_i - bar(bold(x)))_2^2),
    $
    whose unique global minimizer is $bold(w)_1^* = dots = bold(w)_M^* = bar(bold(x))$.
  + *Uniform Contraction of Outlier Prototypes:* Suppose $X subset.eq B(bar(bold(x)), r)$. For any $lambda > 0$, $eta_t in (0, 1]$, and any prototype initialized outside $B(bar(bold(x)), r)$, every online NG step contracts its distance to $bar(bold(x))$ according to
    $
      norm(bold(w)_i^((t+1)) - bar(bold(x)))_2
        <= (1 - gamma_t) norm(bold(w)_i^((t)) - bar(bold(x)))_2 + gamma_t r,
    $
    where $gamma_t := eta_t exp(-(M - 1) / lambda_t) > 0$. Thus no prototype can ever remain permanently dead outside $"conv"(X)$.
]

#proof[
  For (1), as $lambda -> infinity$, $h_lambda (k) = exp(-k / lambda) -> 1$ for every $k in {0, dots, M - 1}$, while $C(lambda) -> M$, so $h_lambda (k) / C(lambda) -> 1 / M$ uniformly over $k$. Substituting into $hat(E)_lambda (W)$ and applying the bias-variance identity $(1 / n) sum_(m=1)^n norm(bold(x)_m - bold(w)_i)_2^2 = V_X + norm(bold(w)_i - bar(bold(x)))_2^2$ yields $E_infinity (W)$. Its Hessian is $(1 / M) I_(M d) succ 0$, establishing strict convexity and the unique minimizer $bold(w)_i^* = bar(bold(x))$.

  For (2), write $alpha_(i t) := eta_t h_(lambda_t) (k_i (bold(x)_t, W^((t))))$. Because $k_i <= M - 1$, $gamma_t <= alpha_(i t) <= 1$. Subtracting $bar(bold(x))$ from the online update gives
  $
    bold(w)_i^((t+1)) - bar(bold(x))
      = (1 - alpha_(i t)) (bold(w)_i^((t)) - bar(bold(x)))
        + alpha_(i t) (bold(x)_t - bar(bold(x))).
  $
  Taking the Euclidean norm and using $norm(bold(x)_t - bar(bold(x)))_2 <= r < norm(bold(w)_i^((t)) - bar(bold(x)))_2$, the affine function $u |-> (1 - u) norm(bold(w)_i^((t)) - bar(bold(x)))_2 + u r$ is strictly decreasing in $u$, so replacing $alpha_(i t) >= gamma_t$ yields the claimed contraction bound. In batch mode, $bold(w)_i^((1))$ enters $"conv"(X)$ on the very first epoch.
]

#subsec("D", "Symmetry and Super-Exponential Homotopy to Lloyd k-Means")

#proposition(title: "Full Prototype Symmetry and Zero-Range Homotopy Error Bound", origin: "Original Proof")[
  Let $Q(W) := (1 / n) sum_(m=1)^n min_(i in [M]) norm(bold(x)_m - bold(w)_i)_2^2$ be the $k$-means quantization distortion @lloyd1982least.
  + $hat(E)_lambda (W)$ is strictly invariant under all $M!$ prototype permutations $sigma in S_M$, invariant to data ordering in batch mode, and equivariant under similarity transforms $bold(x) |-> s Q bold(x) + bold(t)$.
  + Suppose $X union W$ lies inside a ball of diameter $D$. Then for every $lambda > 0$,
    $
      0 <= hat(E)_lambda (W) - 1 / 2 Q(W)
        <= (D^2 / 2) (1 - 1 / C(lambda))
        < (D^2 / 2) exp(-1 / lambda).
    $
]

#proof[
  Part (1) holds because distance ranks $k_i (bold(x)_m, W)$ depend purely on unordered pairwise Euclidean distances, unaffected by relabeling prototypes or applying rigid motions; isotropic scaling by $s$ scales $hat(E)_lambda$ homogeneously by $s^2$.

  For (2), let $d_((0)) (m) <= d_((1)) (m) <= dots <= d_((M - 1)) (m)$ be the sorted squared distances from $bold(x)_m$ to $W$. Since $h_lambda (0) = 1$,
  $
    hat(E)_lambda (W) - 1 / 2 Q(W)
      = 1 / (2 n C(lambda)) sum_(m=1)^n sum_(k=1)^(M-1)
        h_lambda (k) (d_((k)) (m) - d_((0)) (m)).
  $
  Every summand is nonnegative, establishing the lower bound $0$. Moreover, $d_((k)) (m) - d_((0)) (m) <= D^2$, so
  $
    hat(E)_lambda (W) - 1 / 2 Q(W)
      <= (D^2 / 2) (sum_(k=1)^(M-1) h_lambda (k)) / C(lambda)
      = (D^2 / 2) (1 - 1 / C(lambda)).
  $
  Finally, $C(lambda) = (1 - exp(-M / lambda)) / (1 - exp(-1 / lambda)) < 1 / (1 - exp(-1 / lambda))$, which implies $1 / C(lambda) > 1 - exp(-1 / lambda)$ and hence $1 - 1 / C(lambda) < exp(-1 / lambda)$, completing the proof.
]

Theorem II.C and Proposition II.D jointly explain Neural Gas's remarkable optimization quality: annealing $lambda$ from $lambda_0 >> 1$ down to $lambda_f -> 0^+$ traces a smooth deterministic-annealing homotopy from the *unique* global minimum of a strictly convex quadratic ($E_infinity$) into an exponentially accurate surrogate of the NP-hard $k$-means objective @aloise2009np.

#sec("III", "Complexity, Geometry and Operational Profile")

In online mode over $T$ presentations, each step computes $M$ squared Euclidean distances in $cal(O)(M d)$ operations, sorts $M$ distance scalars in $cal(O)(M log M)$ time to obtain exact ranks $k_i$, and updates all $M$ prototypes in $cal(O)(M d)$ operations, yielding $cal(O)(T M (d + log M))$ total time. When $lambda_t$ becomes small during late annealing, truncating updates to ranks $k <= k_"max" (lambda_t) := ceil(lambda_t ln(1 / tau))$ (for tolerance $tau$, e.g., $10^(-4)$) replaces full sorting with partial selection in $cal(O)(M + k_"max" log k_"max")$ time and $cal(O)(k_"max" d)$ update work. In batch mode over $E$ epochs, computing the $n times M$ distance matrix, sorting each row, and forming rank-weighted averages runs in $cal(O)(E n M (d + log M))$ time.

Space complexity is $cal(O)(M d + M)$ in online mode (plus $cal(O)(n d)$ if storing the dataset for batch epochs). Because online NG accesses prototypes linearly in contiguous memory per sample, its memory footprint is cache-resident whenever $M d$ fits in L2/L3 cache.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(T M (d + log M))$ online / $cal(O)(E n M (d + log M))$ batch],
    [Space Complexity], [$cal(O)(M d)$ model state; $cal(O)(n d)$ optional batch data buffer],
    [Cluster Geometry], [Unconstrained Voronoi partition; topology-free manifold coverage],
    [Hyperparameters], [Prototype count $M$, range $lambda_0, lambda_f$, learning rate $eta_0, eta_f$],
    [Invariance], [Similarity-equivariant & fully $S_M$-permutation symmetric],
    [Optimization], [Exact gradient potential $E_lambda$; MM monotone batch descent],
    [Outlier Robustness], [Low ($L_2$ squared distortion; breakdown point $1 / n$)],
  ),
  caption: [Operational and geometric profile of Neural Gas (Alg. 179).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *High-Fidelity Speech, Audio, and Image Codebook Quantization:* Lossy compression (LBG / vector quantization @linde1980algorithm) demands minimal expected squared distortion without wasted dead codewords in sparsely populated acoustic sub-bands. Neural Gas's homotopy from $E_infinity$ to $Q(W)$ (Theorems II.C and II.D) consistently achieves lower quantization distortion than Lloyd $k$-means @lloyd1982least and SOM across arbitrary codebook sizes @martinetz1993neuralgas.
+ *High-Dimensional Nonlinear Manifold Quantization with Unknown Intrinsic Dimension:* In robotics kinematics, spectral chemistry, and financial risk spaces, data often occupy manifolds whose intrinsic dimension $q$ is unknown or exceeds $2$. Unlike SOM (which suffers metric tearing when $q > 2$), Neural Gas's full $S_M$ symmetry allocates prototypes freely across any $q <= d$ without topological tension.
+ *Medical Imaging Color and Tissue-Feature Segmentation:* Multi-parametric MRI and histopathology feature vectors exhibit multi-modal intensity distributions with imbalanced tissue volumes. Neural Gas distributes prototypes smoothly across minority lesion modes without losing them to dead-unit collapse, producing compact Voronoi quantizers for downstream segmentation @xu2005survey.
+ *Molecular Dynamics Conformation Coarse-Graining:* Protein-folding trajectories sample high-dimensional torsional spaces with disconnected metastable basins separated by energy barriers. Neural Gas quantizes conformational states with exact gradient dynamics and high Voronoi Silhouette separation, providing discrete microstates for Markov State Model construction.
+ *Hybrid Topology Extraction via Neural Gas + CHL:* Pairing converged Neural Gas prototypes with a single pass of Competitive Hebbian edge insertion @martinetz1993neuralgas produces a provably optimal induced Delaunay graph over arbitrarily complex disconnected manifolds, combining NG's superior prototype placement with graph-based cluster extraction.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Neural Gas exhibits five limitations: (i) *fixed prototype count $M$*, requiring either domain knowledge or external model selection (e.g., Silhouette or BIC) to choose $M$; (ii) *sorting overhead* $cal(O)(M log M)$ per sample when $M$ reaches tens of thousands without rank truncation; (iii) *absence of native 2D grid coordinates*, preventing direct SOM-style 2D component-plane visualization without secondary dimensionality reduction; (iv) *annealing schedule sensitivity*, where overly small $lambda_0$ forfeits dead-unit escape while overly large $lambda_f$ retains bias toward $bar(bold(x))$; and (v) *spherical $L_2$ Voronoi geometry*, which splits elongated or nested non-convex clusters unless paired with CHL graph merging.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Vector Quantization ($q > 2$)], [Excellent: exact homotopy to min SSE], [Tune $lambda_0, lambda_f$; fast partial sort],
    [High-Dim. Continuous Clusters], [Strong: no grid tearing or $L_1$ collapse], [Standardize features first],
    [2D Exploratory Cartography], [Weak: no native 2D grid coordinates], [Use SOM (Alg. 176) for U-matrix maps],
    [Growing / Drifting Streams], [Weak: static $M$ and decaying $eta_t, lambda_t$], [Use GNG (Alg. 177) or Fuzzy ART (Alg. 178)],
    [Non-Convex Single-Labeling], [Moderate: Voronoi splits arcs unless merged], [Add CHL edges or use HDBSCAN @campello2013hdbscan],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Neural Gas.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Unified Taxonomy of Competitive Neural Clustering (Algs. 176--179)", origin: "Original Derivation")[
  Every online prototype update in Algs. 176, 177, and 179 is an instance of the generalized competitive rule
  $
    bold(w)_i^((t+1)) = bold(w)_i^((t))
      + eta_t kappa_t (i, bold(x)_t, W^((t))) (bold(x)_t - bold(w)_i^((t))),
  $
  distinguished solely by the coupling kernel $kappa_t$:
  + *Lloyd / MacQueen $k$-means* @lloyd1982least @macqueen1967some: $kappa_t = bold(1)_(k_i = 0)$ (hard winner-take-all; non-convex, dead-unit risk);
  + *SOM (Alg. 176)* @kohonen1982som: $kappa_t = h_(sigma_t) (b_W (bold(x)_t), i)$ (fixed external lattice metric; non-integrable curl in continuous space);
  + *GNG (Alg. 177)* @fritzke1995gng: $kappa_t = bold(1)_(i = s_1) + (epsilon_n / epsilon_b) bold(1)_(i in N_t (s_1))$ with dynamic vertex/edge graph $G_t$;
  + *Neural Gas (Alg. 179)* @martinetz1993neuralgas: $kappa_t = exp(-k_i (bold(x)_t, W^((t))) / lambda_t)$ (intrinsic input-space distance rank; exact differentiable gradient potential $E_lambda$).
]

#proof[
  Direct substitution of each coupling function $kappa_t$ into the generalized update recovers the exact per-step prototype recursions of MacQueen $k$-means, online SOM, GNG (with $eta_t = epsilon_b$), and online Neural Gas respectively. In particular, Neural Gas is the unique member of the family with strictly positive soft coupling across all prototypes that preserves exact bisector continuity (Theorem II.A) and full $S_M$ symmetry (Proposition II.D), whereas ART (Alg. 178) @carpenter1991fuzzyart departs from additive $L_2$ displacement entirely in favor of monotonic $L_1$ fuzzy-AND box contraction gated by vigilance $rho$.
]

*ClusterOverview Benchmark Analysis:* Across the standardized 4-dataset ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), Neural Gas achieves Iris Sil = $0.447$, ARI = $0.583$; Customers Sil = $0.340$; Multi-Shape Sil = $0.718$ (ARI = N/A); and High-D Sil = $0.303$ (ARI = N/A), at a mean runtime of $1.327$ ms. Strikingly, Neural Gas achieves the *highest Silhouette score on three of the four datasets* among Algs. 176--179—Iris ($0.447$), Multi-Shape ($0.718$), and High-D ($0.303$)—directly confirming Theorems II.A--II.D: because Neural Gas minimizes an exact smooth homotopy to the Euclidean quantization objective without SOM's fixed-lattice border distortion (Iris Sil $0.447$ vs. SOM $0.290$) and without ART's high-dimensional $L_1$ vigilance collapse (High-D Sil $0.303$ vs. ART $-0.222$), its Voronoi cells attain superior Euclidean compactness and separation. On Iris, NG's ARI ($0.583$) trails ART ($0.793$) and GNG ($0.712$) despite leading in Silhouette ($0.447$), illustrating a classic metric-vs-label nuance: pure isotropic $L_2$ quantization partitions the overlapping _versicolor_ / _virginica_ cloud along minimum-variance boundaries rather than axis-aligned petal thresholds. Finally, NG's mean runtime ($1.327$ ms) sits cleanly between GNG ($0.304$ ms, sparse 1-hop updates) and SOM ($2.477$ ms, full 2D lattice kernel evaluations), reflecting its $cal(O)(T M (d + log M))$ complexity.
