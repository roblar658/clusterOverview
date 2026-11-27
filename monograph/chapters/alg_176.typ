#import "../template/ieee.typ": *

#alg_entry(
  id: 176,
  name: "Self-Organizing Maps (SOM)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "1982",
  paradigm: "Topology-Preserving Competitive Neural Lattice",
  complexity: $cal(O)(E n M d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

The Self-Organizing Map (SOM), introduced by Kohonen @kohonen1982som, couples competitive vector quantization with discrete topological smoothing on a low-dimensional neural lattice. Let $X = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ denote a finite dataset drawn from an underlying probability distribution $P$ on $bb(R)^d$. Fix an undirected lattice graph $G = (V, E_G)$ with $|V| = M$ neurons, embedded in a low-dimensional grid space via fixed coordinates ${bold(r)_i}_(i=1)^M subset bb(R)^q$ (typically a two-dimensional rectangular or hexagonal lattice with $q = 2$, or a one-dimensional chain with $q = 1$), equipped with the lattice metric $d_G (i, j) := norm(bold(r)_i - bold(r)_j)_2$. Each neuron $i in [M]$ carries a prototype vector $bold(w)_i in bb(R)^d$, collected in the codebook matrix $W = (bold(w)_1, dots, bold(w)_M) in bb(R)^(d times M)$.

For any input $bold(x) in bb(R)^d$, Kohonen's Best Matching Unit (BMU) is the nearest prototype in Euclidean metric,
$
  b_W (bold(x)) := arg min_(i in [M]) norm(bold(x) - bold(w)_i)_2^2,
$
with ties broken deterministically by smallest index. The BMU mapping partitions $bb(R)^d$ into $M$ convex Voronoi polyhedra $V_i (W) := {bold(x) in bb(R)^d : b_W (bold(x)) = i}$. Topological coupling across the grid is mediated by a symmetric, nonnegative, unimodal neighborhood kernel $h_sigma : [M] times [M] -> (0, 1]$, canonically the Gaussian lattice kernel
$
  h_sigma (i, j) := exp(- d_G (i, j)^2 / (2 sigma^2)),
$
where $sigma > 0$ is an annealed neighborhood radius satisfying $h_sigma (i, i) = 1$ and $h_sigma (i, j) = h_sigma (j, i)$.

#definition(title: "Topology-Coupled Quantization Distortion Functional")[
  Given dataset $X = {bold(x)_m}_(m=1)^n$, fixed radius $sigma > 0$, codebook $W in bb(R)^(d times M)$, and assignment map $c : [n] -> [M]$, the Heskes-Kohonen coupled distortion functional is
  $
    J_sigma (W, c) := 1 / n sum_(m=1)^n sum_(j=1)^M
      h_sigma (c(m), j) norm(bold(x)_m - bold(w)_j)_2^2.
  $
  Minimizing over $c$ yields the neighborhood-weighted assignment
  $
    c_W^* (m) := arg min_(i in [M]) sum_(j=1)^M
      h_sigma (i, j) norm(bold(x)_m - bold(w)_j)_2^2,
  $
  and the reduced energy $E_sigma (W) := min_c J_sigma (W, c) = J_sigma (W, c_W^*)$.
]

Crucially, Kohonen's classic algorithm uses the Euclidean BMU $b_W (bold(x)_m)$ rather than $c_W^* (m)$. Although both coincide in the zero-neighborhood limit $sigma -> 0^+$ and are close whenever prototypes vary smoothly over the lattice, Section II proves that this seemingly minor distinction separates exact Lyapunov descent (batch Heskes SOM) from non-integrable vector-field dynamics (continuous-space Kohonen SOM).

#algbox(
  title: "Online and Batch Self-Organizing Map (SOM)",
  inputs: [Data $X = {bold(x)_m}_(m=1)^n subset bb(R)^d$; lattice coordinates ${bold(r)_i}_(i=1)^M$; radius schedule $sigma_0 > sigma_f > 0$; step-size schedule $eta_0 > eta_f > 0$ (online) or epochs $E$ (batch).],
  outputs: [Converged prototypes $W = {bold(w)_i}_(i=1)^M$, BMU labels $c(m) = b_W (bold(x)_m)$, and lattice neighborhood map.],
)[
  + *Initialize:* Set $W^((0))$ by linear span along the top-$q$ principal components of $X$ (or random draws from $X$).
  + *Online Mode ($t = 0, dots, T - 1$):*
    - Anneal $sigma_t = sigma_0 (sigma_f / sigma_0)^(t / T)$ and $eta_t = eta_0 (eta_f / eta_0)^(t / T)$.
    - Sample $bold(x)_t in X$ and compute BMU $i^* = arg min_i norm(bold(x)_t - bold(w)_i)_2$.
    - Update every neuron $j in [M]$:
      $
        bold(w)_j <- bold(w)_j + eta_t h_(sigma_t) (i^*, j) (bold(x)_t - bold(w)_j).
      $
  + *Batch Mode ($e = 1, dots, E$):*
    - Anneal $sigma_e$; assign each point $c(m) <- b_W (bold(x)_m)$ (Kohonen) or $c(m) <- c_W^* (m)$ (Heskes).
    - Compute Voronoi cell counts $n_i = |{m : c(m) = i}|$ and sums $bold(s)_i = sum_(c(m) = i) bold(x)_m$ for $i in [M]$.
    - Update each prototype by Nadaraya-Watson lattice smoothing:
      $
        bold(w)_j <- (sum_(i=1)^M h_(sigma_e) (i, j) bold(s)_i) /
          (sum_(i=1)^M h_(sigma_e) (i, j) n_i).
      $
  + *Cluster Extraction:* Return $W$ and labels $c(m) = b_W (bold(x)_m)$, optionally merging lattice cells via U-matrix watershed or hierarchical linkage when $M > k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish four foundational properties of SOM: exact monotonic Lyapunov descent for the Heskes batch formulation, constructive non-integrability of Kohonen's online Euclidean-BMU vector field, absorbing invariance of ordered one-dimensional states, and an exact variance-contraction identity quantifying the SOM border effect.

#subsec("A", "Monotone Descent and Finite Convergence of Batch SOM")

#theorem(title: "Exact Lyapunov Descent of Heskes Batch SOM", origin: "Original Proof")[
  Fix a positive kernel $h_sigma (i, j) > 0$ on $[M]^2$. Alternating between assignment $c^((t+1)) = c_(W^((t)))^*$ and batch prototype smoothing $bold(w)_j^((t+1)) = (sum_(m=1)^n h_sigma (c^((t+1)) (m), j) bold(x)_m) / S_j^((t+1))$, where $S_j^((t+1)) := sum_(m=1)^n h_sigma (c^((t+1)) (m), j) > 0$, produces a monotonically non-increasing energy sequence satisfying
  $
    J_sigma (W^((t)), c^((t))) - J_sigma (W^((t+1)), c^((t+1)))
      >= 1 / n sum_(j=1)^M S_j^((t+1))
        norm(bold(w)_j^((t+1)) - bold(w)_j^((t)))_2^2.
  $
  Moreover, under deterministic tie-breaking, the iteration terminates at a coordinate-wise stationary fixed point in at most $M^n$ steps.
]

#proof[
  Decompose the one-step energy reduction into the assignment phase and the prototype update phase:
  $
    J_sigma (W^((t)), c^((t))) - J_sigma (W^((t+1)), c^((t+1)))
      = Delta_("assign") + Delta_("proto"),
  $
  where $Delta_("assign") := J_sigma (W^((t)), c^((t))) - J_sigma (W^((t)), c^((t+1)))$ and $Delta_("proto") := J_sigma (W^((t)), c^((t+1))) - J_sigma (W^((t+1)), c^((t+1)))$. Because $c^((t+1)) (m)$ minimizes $sum_(j=1)^M h_sigma (i, j) norm(bold(x)_m - bold(w)_j^((t)))_2^2$ independently for each $m in [n]$, we immediately obtain $Delta_("assign") >= 0$.

  For the prototype update step, abbreviate $h_(m j) := h_sigma (c^((t+1)) (m), j)$. For each neuron $j in [M]$, expand the squared residual around $bold(w)_j^((t+1))$:
  $
    norm(bold(x)_m - bold(w)_j^((t)))_2^2
      &= norm(bold(x)_m - bold(w)_j^((t+1)))_2^2
       + norm(bold(w)_j^((t+1)) - bold(w)_j^((t)))_2^2 \
      &quad + 2 chevron.l bold(x)_m - bold(w)_j^((t+1)),
        bold(w)_j^((t+1)) - bold(w)_j^((t)) chevron.r.
  $
  Multiplying by $h_(m j)$ and summing over $m = 1, dots, n$, the weighted cross-term vanishes identically because the batch smoothing rule enforces
  $
    sum_(m=1)^n h_(m j) (bold(x)_m - bold(w)_j^((t+1)))
      = sum_(m=1)^n h_(m j) bold(x)_m - S_j^((t+1)) bold(w)_j^((t+1))
      = bold(0).
  $
  Summing over $j = 1, dots, M$ and dividing by $n$ yields the exact identity
  $
    Delta_("proto") = 1 / n sum_(j=1)^M S_j^((t+1))
      norm(bold(w)_j^((t+1)) - bold(w)_j^((t)))_2^2 >= 0.
  $
  Combining $Delta_("assign") >= 0$ with $Delta_("proto")$ proves the claimed descent inequality. Since $S_j^((t+1)) > 0$, any change $W^((t+1)) != W^((t))$ forces a strict decrease $J_sigma (W^((t+1)), c^((t+1))) < J_sigma (W^((t)), c^((t)))$. Because $W^((t+1))$ is uniquely determined by the discrete assignment $c^((t+1)) in [M]^n$, no assignment can ever repeat after a strict drop. Hence the algorithm must reach $c^((t+1)) = c^((t))$ and $W^((t+1)) = W^((t))$ within at most $|[M]^n| = M^n$ iterations.
]

#subsec("B", "Non-Integrability of Kohonen's Euclidean-BMU Dynamics")

In contrast to the Heskes rule above, Kohonen's standard online SOM with Euclidean BMU selection does not admit a continuous-space Lyapunov energy function.

#theorem(title: "Non-Existence of a Potential for Kohonen's Online SOM", origin: "Original Counterexample")[
  Consider a one-dimensional chain of $M = 3$ prototypes $w_1 < w_2 < w_3$ in $bb(R)$ with symmetric strictly decaying neighborhood weights $h_sigma (i, i) = 1$, $h_1 := h_sigma (1, 2) = h_sigma (2, 3)$, and $h_2 := h_sigma (1, 3)$ satisfying $1 > h_1 > h_2 >= 0$. Let $X ~ P$ have a continuous strictly positive density $p(x) > 0$ on $bb(R)$ with finite first moment. Define the expected Kohonen drift field $F(w_1, w_2, w_3) = (F_1, F_2, F_3)$ by $F_j (W) := bb(E)[h_sigma (b_W (X), j) (X - w_j)]$. Then at every ordered point $w_1 < w_2 < w_3$, the curl component between neurons $1$ and $2$ satisfies
  $
    (partial F_1) / (partial w_2) - (partial F_2) / (partial w_1)
      = ((h_1 - h_2) / 2) (m_(2 3) - w_1) p(m_(2 3)) > 0,
  $
  where $m_(2 3) = (w_2 + w_3) / 2$. Consequently, there exists no $C^2$ potential $U(w_1, w_2, w_3)$ such that $F = -nabla U$.
]

#proof[
  Let $m_(1 2) := (w_1 + w_2) / 2$ and $m_(2 3) := (w_2 + w_3) / 2$ be the two Voronoi boundaries separating $V_1 = (-infinity, m_(1 2))$, $V_2 = (m_(1 2), m_(2 3))$, and $V_3 = (m_(2 3), infinity)$. Writing out the first two drift components explicitly gives
  $
    F_1 &= integral_(-infinity)^(m_(1 2)) (x - w_1) p(x) dif x
      + h_1 integral_(m_(1 2))^(m_(2 3)) (x - w_1) p(x) dif x \
      &quad + h_2 integral_(m_(2 3))^infinity (x - w_1) p(x) dif x, \
    F_2 &= h_1 integral_(-infinity)^(m_(1 2)) (x - w_2) p(x) dif x
      + integral_(m_(1 2))^(m_(2 3)) (x - w_2) p(x) dif x \
      &quad + h_1 integral_(m_(2 3))^infinity (x - w_2) p(x) dif x.
  $
  Because $w_2$ enters $F_1$ solely through the limits $m_(1 2)$ and $m_(2 3)$ with $partial m_(1 2) / partial w_2 = partial m_(2 3) / partial w_2 = 1 / 2$, the Leibniz integral rule yields
  $
    (partial F_1) / (partial w_2)
      &= ((1 - h_1) / 2) (m_(1 2) - w_1) p(m_(1 2)) \
      &quad + ((h_1 - h_2) / 2) (m_(2 3) - w_1) p(m_(2 3)).
  $
  Symmetrically, $w_1$ enters $F_2$ solely through $m_(1 2)$ with $partial m_(1 2) / partial w_1 = 1 / 2$, giving
  $
    (partial F_2) / (partial w_1)
      = ((h_1 - 1) / 2) (m_(1 2) - w_2) p(m_(1 2)).
  $
  Observe that $m_(1 2) - w_1 = (w_2 - w_1) / 2 = -(m_(1 2) - w_2)$, so the $p(m_(1 2))$ terms in $(partial F_1) / (partial w_2)$ and $(partial F_2) / (partial w_1)$ are algebraically identical and cancel upon subtraction:
  $
    (partial F_1) / (partial w_2) - (partial F_2) / (partial w_1)
      = ((h_1 - h_2) / 2) (m_(2 3) - w_1) p(m_(2 3)).
  $
  Since $w_1 < w_2 < w_3$, we have $m_(2 3) > w_2 > w_1$, while $h_1 > h_2$ and $p(m_(2 3)) > 0$ by assumption; thus the curl is strictly positive everywhere on the ordered cone. If $F = -nabla U$ held for some $C^2$ potential $U$, Schwarz's theorem would require $(partial F_1) / (partial w_2) = -(partial^2 U) / (partial w_2 partial w_1) = (partial F_2) / (partial w_1)$, a contradiction. Notice that as $sigma -> 0^+$, $h_1 = h_2 = 0$, restoring zero curl and recovering the exact $k$-means potential.
]

#subsec("C", "Absorbing Invariance of Ordered One-Dimensional States")

Despite the absence of a global potential, Kohonen's online update possesses strong topological self-organization guarantees in one dimension @kohonen1982som.

#theorem(title: "Positive Invariance of Ordered Configurations in 1D SOM", origin: "Classical — Cited")[
  Let $G$ be a one-dimensional chain of $M$ neurons with monotone non-increasing neighborhood weights $1 = h_0 >= h_1 >= dots >= h_(M-1) >= 0$, where $h_(|i - j|) := h_sigma (i, j)$, and step size $eta in (0, 1)$. Then the strictly ordered simplex $Omega_+ := {W in bb(R)^M : w_1 < w_2 < dots < w_M}$ (and symmetrically $Omega_-$) is positively invariant under every online Kohonen update @kohonen1982som: if $W in Omega_+$, then $W' in Omega_+$ for every input $x in bb(R)$.
]

#proof[
  Let $W in Omega_+$, draw any $x in bb(R)$, and let $i^* = b_W (x)$ be the Euclidean BMU. Fix any adjacent index pair $(j, j + 1)$ with $1 <= j <= M - 1$, and denote $a := eta h_(|i^* - j|)$ and $b := eta h_(|i^* - (j + 1)|)$. Since $eta in (0, 1)$ and $h_r in [0, 1]$, both $a$ and $b$ lie in $[0, 1)$. The updated gap is
  $
    w'_(j+1) - w'_j
      = (1 - b) w_(j+1) - (1 - a) w_j + (b - a) x.
  $
  We examine the two exhaustive cases for the position of $i^*$:
  - *Case 1 ($j + 1 <= i^*$):* Then $|i^* - (j + 1)| < |i^* - j|$, so monotonicity of $h_r$ implies $0 <= a <= b < 1$. Rewriting the gap as
    $
      w'_(j+1) - w'_j = (1 - b)(w_(j+1) - w_j) + (b - a)(x - w_j),
    $
    we note that $w_j < w_(j+1) <= w_(i^*)$. Because $i^*$ is the BMU, $|x - w_(i^*)| <= |x - w_j|$, which implies $x >= (w_j + w_(i^*)) / 2 > w_j$. Thus $(1 - b)(w_(j+1) - w_j) > 0$ and $(b - a)(x - w_j) >= 0$, yielding $w'_(j+1) > w'_j$.
  - *Case 2 ($j >= i^*$):* Then $|i^* - j| < |i^* - (j + 1)|$, so $0 <= b <= a < 1$. Rewriting as
    $
      w'_(j+1) - w'_j = (1 - a)(w_(j+1) - w_j) + (a - b)(w_(j+1) - x),
    $
    BMU optimality $|x - w_(i^*)| <= |x - w_(j+1)|$ with $w_(i^*) <= w_j < w_(j+1)$ implies $x <= (w_(i^*) + w_(j+1)) / 2 < w_(j+1)$. Hence $(1 - a)(w_(j+1) - w_j) > 0$ and $(a - b)(w_(j+1) - x) >= 0$, again giving $w'_(j+1) > w'_j$.

  Since this holds for all $j in {1, dots, M - 1}$, $W' in Omega_+$.
]

#subsec("D", "Convex Hull Containment and Exact Border-Contraction Identity")

#proposition(title: "Equivariance, Convex Hull Invariance, and Exact Variance Contraction", origin: "Original Proof")[
  For both online and batch SOM with $eta_t in [0, 1]$ and $h_sigma (i, j) in [0, 1]$:
  + *Similarity Equivariance:* Transforming $bold(x)_m |-> s Q bold(x)_m + bold(t)$ ($s > 0$, $Q in O(d)$, $bold(t) in bb(R)^d$) transforms prototypes identically as $bold(w)_j |-> s Q bold(w)_j + bold(t)$, whereas anisotropic feature scaling breaks equivariance.
  + *Convex Hull Invariance:* If $bold(w)_j^((0)) in "conv"(X)$ for all $j$, then $bold(w)_j^((t)) in "conv"(X)$ for all $t >= 0$.
  + *Exact Variance Contraction (Border Effect):* At any batch fixed point $W^*$, let $bar(bold(x))_i$ and $n_i$ be the empirical centroid and count of Voronoi cell $i$, let $bar(bold(x)) = (1 / n) sum_m bold(x)_m$, let $S_j := sum_(i=1)^M n_i h_sigma (i, j)$, and let $H_i := sum_(j=1)^M h_sigma (i, j)$. Then
    $
      sum_(j=1)^M S_j norm(bold(w)_j^* - bar(bold(x)))_2^2
        &= sum_(i=1)^M n_i H_i norm(bar(bold(x))_i - bar(bold(x)))_2^2 \
        &quad - sum_(j=1)^M sum_(i=1)^M n_i h_sigma (i, j)
          norm(bar(bold(x))_i - bold(w)_j^*)_2^2.
    $
]

#proof[
  Property (1) follows directly from preservation of Euclidean distance ratios under similarity transforms, leaving BMU indices invariant, and linearity of the convex-combination prototype updates. For (2), each online update $bold(w)'_j = (1 - eta_t h) bold(w)_j + (eta_t h) bold(x)_t$ has coefficient $eta_t h in [0, 1]$, and each batch update is a normalized nonnegative linear combination of $bold(x)_m in X$; by induction both preserve $"conv"(X)$.

  For (3), at a batch fixed point $bold(w)_j^* = (1 / S_j) sum_(i=1)^M n_i h_sigma (i, j) bar(bold(x))_i$. For each fixed $j$, expand $norm(bar(bold(x))_i - bar(bold(x)))_2^2$ around $bold(w)_j^*$:
  $
    norm(bar(bold(x))_i - bar(bold(x)))_2^2
      &= norm(bar(bold(x))_i - bold(w)_j^*)_2^2
       + norm(bold(w)_j^* - bar(bold(x)))_2^2 \
      &quad + 2 chevron.l bar(bold(x))_i - bold(w)_j^*,
        bold(w)_j^* - bar(bold(x)) chevron.r.
  $
  Weighting by $n_i h_sigma (i, j)$ and summing over $i = 1, dots, M$ annihilates the cross-term because $sum_i n_i h_sigma (i, j) (bar(bold(x))_i - bold(w)_j^*) = bold(0)$. Summing over $j = 1, dots, M$ and rearranging produces the exact identity. On a vertex-transitive (toroidal) lattice where $H_i equiv H$, dividing by $H$ shows that the effective between-prototype dispersion is strictly smaller than the Voronoi-centroid dispersion whenever neighboring cells have distinct centroids; on a bounded planar grid, interior pulls contract border prototypes strictly inward toward $"conv"(X)$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

In online mode over $T$ sample presentations, each step evaluates $M$ squared Euclidean distances in $cal(O)(M d)$ operations and updates up to $M$ prototypes in $cal(O)(M d)$ operations, yielding $cal(O)(T M d)$ total time (or $cal(O)(E n M d)$ for $E$ full passes over $n$ points). In batch mode, an essential two-stage factorization avoids $cal(O)(n M d)$ smoothing costs: after computing BMU assignments in $cal(O)(n M d)$ and accumulating cell sums $bold(s)_i$ in $cal(O)(n d)$, multiplying the $M times M$ kernel matrix $H_sigma$ against the $M times d$ centroid matrix costs only $cal(O)(M^2 d)$ per epoch. Total batch complexity is therefore $cal(O)(E (n M d + M^2 d))$, dominated by $cal(O)(E n M d)$ whenever $M <= n$. Space complexity is $cal(O)(M d + M^2)$ for prototypes and the precomputed grid kernel, plus $cal(O)(n d)$ if the dataset resides in memory.

Because prototype matrix-vector products are dense BLAS Level-3 GEMM operations in batch mode (computing pairwise distance matrices $X W$ in blocked tiles), batch SOM achieves substantially higher cache locality and SIMD utilization than single-sample online updates, though online SOM supports true streaming with $cal(O)(M d)$ memory.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(E n M d)$ batch / $cal(O)(T M d)$ online; GEMM-friendly in batch mode],
    [Space Complexity], [$cal(O)(M d + M^2)$ model state; $cal(O)(n d)$ optional batch buffer],
    [Cluster Geometry], [Convex Voronoi cells linked by fixed $q$-D lattice topology],
    [Hyperparameters], [Lattice shape $(M_1, M_2)$, topology (rect / hex), $sigma_0, sigma_f$, $eta_0, eta_f$],
    [Invariance], [Rigid-motion & isotropic-scale equivariant; feature-scale sensitive],
    [Determinism], [Deterministic under PCA init + batch mode; order-sensitive online],
    [Outlier Robustness], [Low (squared $L_2$ pull and mean smoothing; breakdown point $1 / n$)],
  ),
  caption: [Operational and geometric profile of Self-Organizing Maps (Alg. 176).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *High-Dimensional Transcriptomics and Gene-Expression Cartography:* When profiling thousands of genes across cell states, biologists need simultaneous clustering and an interpretable 2D visual atlas. SOM projects continuous gene-expression gradients onto adjacent hexagonal map units, preserving local transcriptional neighborhoods while summarizing expression profiles via component planes @xu2005survey.
+ *Industrial Process Monitoring and Fault-Trajectory Visualization:* Multivariate chemical reactors and power turbines generate correlated sensor streams moving smoothly through operating regimes. Mapping trajectories onto a trained SOM with its unified distance matrix (U-matrix) exposes regime boundaries as high-distance ridges and flags incipient drift as movement across low-density grid cells.
+ *Hyperspectral Remote Sensing and Land-Cover Quantization:* Airborne imaging spectrometers yield $d in [100, 400]$ spectral bands per pixel with continuous physical mixing between endmembers. SOM's smooth lattice interpolation models spectral mixing transitions more faithfully than unstructured vector quantizers while compressing scenes into compact spectral codebooks @jain2010data.
+ *Macroeconomic Regime and Yield-Curve Topology:* Sovereign yield curves and macro-financial indicators evolve along low-dimensional nonlinear manifolds (level, slope, curvature). A 2D SOM preserves ordinal proximity between crisis, expansion, and tightening regimes, enabling interpretable state-transition heatmaps across historical cycles.
+ *Biomedical EEG / ECG Biosignal State Mapping:* Continuous physiological monitoring produces noisy multi-channel spectral features whose pathological states form connected submanifolds. Topological smoothing regularizes prototype estimation across sparsely populated clinical sub-states, preventing isolated overfitting while exposing patient-specific state trajectories.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SOM exhibits four systematic mathematical failure modes: (i) *topological defects and folding* ("butterfly twists") when the intrinsic data manifold dimension exceeds the lattice dimension $q$ or when $sigma_t$ decays too rapidly before global unfolding; (ii) *boundary shrinkage* quantified by Proposition II.D, which biases edge prototypes inward and inflates quantization error near support boundaries; (iii) *fixed-topology mismatch*, where disconnected clusters force intermediate "interpolating" neurons into empty regions between modes (often lowering Silhouette scores); and (iv) *feature-scale dominance*, where unstandardized coordinates distort Euclidean BMU selection.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [2D Visual EDA & Atlases], [Strong: native grid layout & U-matrix], [Border shrinkage; toroidal grid],
    [Continuous Manifolds ($q <= 2$)], [Strong: smooth topological sheet], [Fold risk; use PCA init],
    [Disconnected / Non-Convex], [Weak: dead bridge nodes across gaps], [Use GNG (Alg. 177) or HDBSCAN @campello2013hdbscan],
    [Unknown Intrinsic Dim. ($q > 2$)], [Poor: metric tearing on 2D grid], [Use Neural Gas (Alg. 179)],
    [Streaming Concept Drift], [Moderate: fixed $M$, forgetting if $eta$ high], [Use Fuzzy ART (Alg. 178) or DenStream @cao2006denstream],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for SOM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

#proposition(title: "Zero-Neighborhood Reduction of SOM to Lloyd and MacQueen k-Means", origin: "Original Derivation")[
  Let $d_(min) := min_(i != j) d_G (i, j) > 0$ be the minimum inter-neuron lattice spacing. As $sigma -> 0^+$, $h_sigma (i, j) -> bold(1)_(i = j)$ at rate $cal(O)(exp(-d_(min)^2 / (2 sigma^2)))$. Consequently:
  + $J_sigma (W, c)$ converges uniformly on bounded prototype sets to the $k$-means objective $(1 / n) sum_(m=1)^n norm(bold(x)_m - bold(w)_(c(m)))_2^2$ @lloyd1982least;
  + Batch SOM (both Kohonen and Heskes rules) reduces identically to Lloyd's alternating $k$-means algorithm @lloyd1982least;
  + Online SOM reduces identically to MacQueen's stochastic online $k$-means @macqueen1967some @bottou1995convergence.
]

#proof[
  For any $i != j$, $0 < h_sigma (i, j) <= exp(-d_(min)^2 / (2 sigma^2)) -> 0$ as $sigma -> 0^+$, whereas $h_sigma (i, i) = 1$. Substituting $h_sigma (i, j) = bold(1)_(i = j)$ into $J_sigma (W, c)$ collapses the inner sum to the single term $j = c(m)$, recovering the exact $k$-means mean squared error. Similarly, both $c_W^* (m)$ and $b_W (bold(x)_m)$ collapse to $arg min_i norm(bold(x)_m - bold(w)_i)_2^2$, the batch numerator/denominator reduce to the Voronoi cell mean $bold(s)_j / n_j$, and the online update modifies solely the BMU via $bold(w)_(i^*) <- bold(w)_(i^*) + eta_t (bold(x)_t - bold(w)_(i^*))$.
]

Within Category VIII, SOM serves as the classical fixed-lattice prototype baseline @xu2005survey @hastie2009elements. Removing the fixed grid metric $d_G (i, j)$ in favor of input-space distance ranks yields Neural Gas (Alg. 179) @martinetz1993neuralgas, which restores an exact differentiable energy potential; allowing dynamic node insertion and Competitive Hebbian edge creation yields Growing Neural Gas (Alg. 177) @fritzke1995gng; adding vigilance-gated category creation yields Adaptive Resonance Theory (Alg. 178) @carpenter1987art.

*ClusterOverview Benchmark Analysis:* Across the standardized 4-dataset ClusterOverview suite (evaluated with Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing), SOM achieves Iris Sil = $0.290$, ARI = $0.562$; Customers Sil = $0.315$; Multi-Shape Sil = $0.712$ (ARI = N/A); and High-D Sil = $0.301$ (ARI = N/A), at a mean runtime of $2.477$ ms. These measurements align precisely with our theoretical derivations. On Multi-Shape, SOM attains a strong Silhouette of $0.712$ because its smooth lattice aligns cleanly with well-separated geometric lobes. Conversely, on Iris, SOM records the lowest Silhouette ($0.290$) and ARI ($0.562$) among Algs. 176--179: Proposition II.D explains how positive neighborhood coupling pulls adjacent grid prototypes across the overlapping _Iris versicolor_ / _virginica_ decision boundary, inducing border contraction and mixed-species cells. On High-D, SOM maintains a solid positive Silhouette ($0.301$), avoiding the high-dimensional hyperrectangle overlap that collapses Fuzzy ART ($-0.222$), while its $2.477$ ms mean runtime is the highest of the four neural algorithms because every sample presentation evaluates exponential neighborhood weights and updates all $M$ grid prototypes.
