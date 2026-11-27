#import "../template/ieee.typ": *

#alg_entry(
  id: 152,
  name: "Possibilistic C-Means (PCM)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "1993",
  paradigm: "Decoupled Typicality Mode-Seeking",
  complexity: $O(n c d T)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Fuzzy C-Means (Alg. 151 @bezdek1981fcm) enforces the probabilistic sum-to-one constraint $sum_(j=1)^c u_(i j) = 1$, forcing even extreme outliers equidistant from all prototypes to receive substantial membership $1 / c$ and corrupt every cluster center. Krishnapuram and Keller @krishnapuram1993pcm resolved this structural vulnerability by interpreting cluster membership as _possibility_ (or _typicality_) rather than relative probability: a point's degree of belonging $t_(i j) in [0, 1]$ to cluster $j$ should depend solely on its absolute compatibility with prototype $bold(v)_j$, independent of all other clusters.

Simply dropping the simplex constraint from FCM yields the trivial zero solution $t_(i j) = 0$. PCM therefore augments the distortion functional with a complementary penalty that rewards high typicalities while balancing cluster-specific scale parameters $bold(eta) = (eta_1, dots, eta_c)^top in RR_(>0)^c$. Given $X = {bold(x)_1, dots, bold(x)_n} subset RR^d$, an SPD metric matrix $A succ 0$ with $d_(i j)^2 = norm(bold(x)_i - bold(v)_j)_A^2$, and fuzzifier $m > 1$, PCM minimizes
$
  J_(m, bold(eta)) (T, V) &= sum_(j=1)^c sum_(i=1)^n t_(i j)^m d_(i j)^2 \
                          &+ sum_(j=1)^c eta_j sum_(i=1)^n (1 - t_(i j))^m,
$
over typicality matrices $T = [t_(i j)] in [0, 1]^(n times c)$ and prototypes $V = (bold(v)_1, dots, bold(v)_c) in RR^(d times c)$.

#definition(title: "Possibilistic c-Partition Space")[
  The _possibilistic $c$-partition space_ $M_(p c) subset [0, 1]^(n times c)$ consists of all matrices $T = [t_(i j)]$ satisfying:
  1. _Interval typicality_: $t_(i j) in [0, 1]$ for all $i in {1, dots, n}$ and $j in {1, dots, c}$;
  2. _Non-empty clusters_: $0 < sum_(i=1)^n t_(i j) < n$ for every cluster $j$;
  3. _Pointwise support_: $max_(1 <= j <= c) t_(i j) > 0$ for every observation $i$ (satisfied automatically for finite $d_(i j)$).
]

The bandwidth parameter $eta_j > 0$ sets the squared distance at which typicality drops to $1 / 2$. Because $bold(eta)$ cannot be optimized freely inside $J_(m, bold(eta))$ without collapsing to $0$, PCM initializes $V^((0))$ via FCM and estimates $eta_j$ from the initial fuzzy intra-cluster scatter:
$
  eta_j = K dot (sum_(i=1)^n (u_(i j)^((0)))^m norm(bold(x)_i - bold(v)_j^((0)))_A^2) / (sum_(i=1)^n (u_(i j)^((0)))^m),
$
typically with $K = 1$, optionally refined after a first pass via $alpha$-cut averages @krishnapuram1993pcm.

#algbox(
  title: "Possibilistic C-Means (PCM)",
  inputs: [Data $X = {bold(x)_i}_(i=1)^n subset RR^d$, clusters $c$, fuzzifier $m > 1$, metric $A succ 0$, factor $K > 0$, tolerance $epsilon > 0$, max iterations $T_("max")$.],
  outputs: [Typicality matrix $T^* in M_(p c)$ and prototypes $V^* = {bold(v)_j^*}_(j=1)^c$.],
)[
  1. *Warm-start:* Run FCM (Alg. 151) to obtain initial $(U^((0)), V^((0)))$.
  2. *Estimate scales:* For $j = 1, dots, c$, compute $eta_j$ from $(U^((0)), V^((0)))$ with scaling factor $K$ (guard $eta_j >= eta_("min") > 0$).
  3. Set $t <- 0$ and compute $T^((0))$ from $V^((0))$ via Step 5.
  4. *Repeat* until $norm(V^((t)) - V^((t-1)))_F < epsilon$ or $t = T_("max")$:
  5. #h(1.0em) *Typicality update:* For $i = 1, dots, n$ and $j = 1, dots, c$,
     $ t_(i j)^((t+1)) = (1 + (d_(i j)^2 / eta_j)^(1 / (m - 1)))^(-1), $
     where $d_(i j)^2 = norm(bold(x)_i - bold(v)_j^((t)))_A^2$.
  6. #h(1.0em) *Prototype update:* For $j = 1, dots, c$,
     $ bold(v)_j^((t+1)) = (sum_(i=1)^n (t_(i j)^((t+1)))^m bold(x)_i) / (sum_(i=1)^n (t_(i j)^((t+1)))^m). $
  7. #h(1.0em) Set $t <- t + 1$.
  8. *Return* $(T^((t)), V^((t)))$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

PCM's strengths (redescending outlier resistance) and principal pathology (coincident cluster collapse) both emerge directly from the exact mathematical decoupling of its objective across clusters.

#subsec("A", "Elementwise Separability and Monotone Descent")

#lemma(title: "Exact Decoupled Typicality and AO Descent", origin: "Original Derivation")[
  Let $m > 1$, $A succ 0$, and $bold(eta) in RR_(>0)^c$.
  1. For fixed $V$, $T mapsto J_(m, bold(eta)) (T, V)$ splits into $n c$ independent strictly convex scalar problems on $[0, 1]$, each having the unique global minimizer
     $ t_(i j)^* = (1 + (d_(i j)^2 / eta_j)^(1 / (m - 1)))^(-1), $
     with $t_(i j)^* = 1$ when $d_(i j) = 0$ and $t_(i j)^* = 1 / 2$ precisely on the ellipsoid $d_(i j)^2 = eta_j$.
  2. For fixed $T$ with $sum_(i=1)^n t_(i j)^m > 0$, $V mapsto J_(m, bold(eta)) (T, V)$ is strictly convex quadratic with unique minimizer $bold(v)_j^* = (sum_i t_(i j)^m bold(x)_i) / (sum_i t_(i j)^m)$.
  3. Alternating optimization monotonically decreases $J_(m, bold(eta)) (T^((t)), V^((t)))$ and converges subsequentially to KKT stationary points in $[0, 1]^(n times c) times "conv"(X)^c$.
]

#proof[
  Fix $(i, j)$ and let $g_(i j) (t) = t^m d_(i j)^2 + eta_j (1 - t)^m$ for $t in [0, 1]$. Its derivatives on $(0, 1)$ are
  $
    g'_(i j) (t) &= m t^(m - 1) d_(i j)^2 - m eta_j (1 - t)^(m - 1), \
    g''_(i j) (t) &= m (m - 1) [t^(m - 2) d_(i j)^2 + eta_j (1 - t)^(m - 2)].
  $
  Since $m > 1$, $d_(i j)^2 >= 0$, and $eta_j > 0$, we have $g''_(i j) (t) > 0$ strictly on $(0, 1)$, so $g_(i j)$ is strictly convex on $[0, 1]$. If $d_(i j) = 0$, then $g_(i j) (t) = eta_j (1 - t)^m$ is uniquely minimized at $t_(i j)^* = 1$. If $d_(i j) > 0$, then $g'_(i j) (0^+) = -m eta_j < 0$ and $g'_(i j) (1^-) = m d_(i j)^2 > 0$, guaranteeing a unique interior root solving
  $ ((1 - t) / t)^(m - 1) = d_(i j)^2 / eta_j => t_(i j)^* = (1 + (d_(i j)^2 / eta_j)^(1 / (m - 1)))^(-1). $
  Setting $d_(i j)^2 = eta_j$ gives $t_(i j)^* = (1 + 1)^(-1) = 1 / 2$. Part (2) follows identically to Lemma 1 of FCM since $t_(i j)^* > 0$ for all finite distances ensures $sum_(i=1)^n (t_(i j)^*)^m > 0$. Combining exact strict block minimizations yields $J_(m, bold(eta)) (T^((t+1)), V^((t+1))) <= J_(m, bold(eta)) (T^((t)), V^((t)))$ with strict descent off fixed points, bounded below by $0$ inside the compact product $[0, 1]^(n times c) times "conv"(X)^c$.
]

#subsec("B", "Redescending M-Estimator Equivalence and Exact Influence Bound")

#theorem(title: "Reduction to a Strictly Redescending M-Estimator", origin: "Original Proof")[
  Eliminating $T$ from $J_(m, bold(eta))$ yields the exact reduced prototype objective
  $
    R_(m, bold(eta)) (V) = min_(T in [0, 1]^(n times c)) J_(m, bold(eta)) (T, V) = sum_(j=1)^c sum_(i=1)^n rho_(m, eta_j) (d_(i j)^2),
  $
  where the induced robust loss is
  $
    rho_(m, eta) (r) = (r^(-1 / (m - 1)) + eta^(-1 / (m - 1)))^(1 - m) = r (1 + (r / eta)^(1 / (m - 1)))^(1 - m).
  $
  In particular:
  1. For $m = 2$, $rho_(2, eta) (r) = (eta r) / (eta + r)$ is exactly the classical _Geman-McClure_ robust M-estimator loss;
  2. For every $m > 1$, $rho_(m, eta)$ is bounded above by $eta$ ($lim_(r -> infinity) rho_(m, eta) (r) = eta$) and agrees locally with squared error ($rho_(m, eta) (r) = r + o(r)$ as $r -> 0^+$);
  3. For $A = I_d$, the Euclidean influence vector $psi(bold(x); bold(v)) = -nabla_(bold(v)) rho_(m, eta) (norm(bold(x) - bold(v))_2^2)$ is strictly redescending, vanishing at rate $norm(psi)_2 = O(d^(-(m + 1) / (m - 1)))$ as $d = norm(bold(x) - bold(v))_2 -> infinity$, with exact finite gross-error sensitivity
     $ sup_(d >= 0) norm(psi)_2 = 2 sqrt(eta) ((m - 1) / (m + 1))^((m - 1) / 2) ((m + 1) / (2 m))^m. $
]

#proof[
  Let $q = 1 / (m - 1) > 0$ and $w = (r / eta)^q$, so $r = eta w^(m - 1)$. By Lemma 1, $t^* = (1 + w)^(-1)$ and $1 - t^* = w (1 + w)^(-1)$. Substituting into $rho_(m, eta) (r) = (t^*)^m r + eta (1 - t^*)^m$ produces
  $
    rho_(m, eta) (r) &= (1 + w)^(-m) eta w^(m - 1) + eta w^m (1 + w)^(-m) \
                     &= eta w^(m - 1) (1 + w)^(1 - m) = r (1 + (r / eta)^q)^(1 - m).
  $
  Factoring $r = (r^(-q))^(1 - m)$ gives $(r^(-q) + eta^(-q))^(1 - m)$. At $m = 2$ ($q = 1$), this simplifies immediately to $(r^(-1) + eta^(-1))^(-1) = (eta r) / (eta + r)$. Writing $rho_(m, eta) (r) = eta (1 + (eta / r)^q)^(1 - m)$ verifies $lim_(r -> infinity) rho_(m, eta) (r) = eta$, while $r (1 + (r / eta)^q)^(1 - m) = r + o(r)$ as $r -> 0^+$.

  By the envelope theorem, $rho'_(m, eta) (r) = (t^*(r))^m = (1 + (r / eta)^q)^(-m)$. Hence for $A = I_d$ and $d = norm(bold(x) - bold(v))_2$,
  $ psi(bold(x); bold(v)) = 2 (1 + (d^2 / eta)^q)^(-m) (bold(x) - bold(v)). $
  Its magnitude is $h(d) = 2 d (1 + (d^2 / eta)^q)^(-m)$. As $d -> infinity$, $h(d) tilde.op 2 eta^(m q) d^(1 - 2 m q) = 2 eta^(m / (m - 1)) d^(-(m + 1) / (m - 1)) -> 0$ because $(m + 1) / (m - 1) > 1$ for every $m > 1$. To compute $sup_(d >= 0) h(d)$, differentiate $ln h(d)$ with respect to $d$:
  $ (dif) / (dif d) ln h(d) = 1 / d - (2 m q (d^2 / eta)^q) / (d (1 + (d^2 / eta)^q)) = 0. $
  Setting $w = (d^2 / eta)^q$ gives $1 + w = 2 m q w$, hence $w^* = 1 / (2 m q - 1) = (m - 1) / (m + 1)$, which uniquely maximizes $h(d)$ since $h(0) = h(infinity) = 0$. Substituting $d^* = sqrt(eta) (w^*)^((m - 1) / 2)$ and $1 + w^* = 2 m / (m + 1)$ yields the exact gross-error bound.
]

#subsec("C", "Global Coincidence Degeneracy")

#theorem(title: "Coincident-Cluster Collapse of Global Optima", origin: "Original Proof")[
  Define the single-cluster mode functional $F_(m, eta) (bold(v)) = sum_(i=1)^n rho_(m, eta) (norm(bold(x)_i - bold(v))_A^2)$. Then the reduced PCM objective decouples completely as
  $ R_(m, bold(eta)) (V) = sum_(j=1)^c F_(m, eta_j) (bold(v)_j), $
  with zero cross-cluster interaction. Consequently:
  1. Every global minimizer $V^*$ satisfies $bold(v)_j^* in arg min_(bold(v) in RR^d) F_(m, eta_j) (bold(v))$ independently for each $j$;
  2. Whenever $eta_1 = dots = eta_c = eta$ and $F_(m, eta)$ attains a unique global minimum at $bold(v)^*$, the *unique global minimizer* of PCM is completely degenerate: $bold(v)_1^* = bold(v)_2^* = dots = bold(v)_c^* = bold(v)^*$.
]

#proof[
  The additive decomposition $R_(m, bold(eta)) (V) = sum_(j=1)^c F_(m, eta_j) (bold(v)_j)$ follows immediately from Theorem 1 because each summand depends exclusively on the single vector $bold(v)_j$ and no constraint couples $bold(v)_j$ to $bold(v)_l$ ($j != l$). Therefore,
  $ min_(V in RR^(d times c)) R_(m, bold(eta)) (V) = sum_(j=1)^c min_(bold(v)_j in RR^d) F_(m, eta_j) (bold(v)_j), $
  establishing (1). When $eta_1 = dots = eta_c = eta$, all $c$ subproblems are identical; uniqueness of $bold(v)^* = arg min F_(m, eta)$ immediately forces $bold(v)_j^* = bold(v)^*$ for every $j = 1, dots, c$, proving (2).
]

Theorem 2 formalizes why PCM operates fundamentally as $c$ parallel _mode-seekers_ @krishnapuram1993pcm rather than a partitioner: distinct clusters survive only if each prototype remains trapped inside its own local basin of $F_(m, eta_j)$. If two clusters overlap or if $eta_j$ is slightly oversized, a smaller cluster's local basin vanishes and its prototype slides onto a neighboring higher-density mode.

#subsec("D", "Outlier Rejection and Local Breakdown Stability")

#proposition(title: "Resolution of the Equidistant Noise Paradox and Local Stability", origin: "Original Proof")[
  Consider $n_0 >= 2$ inliers inside a ball $B(bold(0), r)$ with $eta > 2 r^2$ and $A = I_d$, augmented by $k$ arbitrary outliers ${bold(y)_s}_(s=1)^k$ at distance $norm(bold(y)_s)_2 >= D$.
  1. In the equidistant noise setup of Alg. 151 ($bold(x)_0 = (0, L)^top$, prototypes $(plus.minus a, 0)^top$), PCM assigns $t_(0 j) = O(L^(-2 / (m - 1))) -> 0$ and induces a vanishing prototype shift $O(L^(-(m + 1) / (m - 1))) -> 0$ as $L -> infinity$.
  2. For *any* outlier count $k$ (even $k >> n_0$, i.e., contamination fraction $k / (n_0 + k) -> 1$), there exists a finite threshold $D_0 (k, n_0, eta, m, r)$ such that whenever $D >= D_0$, $F_(m, eta)$ retains a strict local minimizer inside $B(bold(0), 2 r)$ regardless of outlier locations.
]

#proof[
  Part (1) follows directly from Lemma 1 ($t_(0 j) = (1 + ((a^2 + L^2) / eta_j)^(1/(m-1)))^(-1) = O(L^(-2/(m-1)))$) and Theorem 1(3), which bounds the shift numerator by $t_(0 j)^m L = O(L^(-(m+1)/(m-1)))$.

  For part (2), examine the outward radial directional derivative of $F_(m, eta) (bold(v))$ at any boundary point $bold(v) in partial B(bold(0), 2 r)$, with unit outward normal $bold(u) = bold(v) / (2 r)$. For every inlier $bold(x)_i in B(bold(0), r)$, we have $chevron.l bold(v) - bold(x)_i, bold(u) chevron.r >= 2 r - r = r$ and $norm(bold(v) - bold(x)_i)_2 <= 3 r$, so each inlier's typicality satisfies $t_(i) >= t_("min") = (1 + (9 r^2 / eta)^q)^(-1) > 0$. Thus the aggregate inward gradient from the $n_0$ inliers satisfies
  $ chevron.l nabla F_("in") (bold(v)), bold(u) chevron.r >= 2 n_0 t_("min")^m r > 0. $
  Meanwhile, every outlier satisfies $norm(bold(v) - bold(y)_s)_2 >= D - 2 r$, so by Theorem 1(3) the total outlier pull is bounded in magnitude by
  $ norm(nabla F_("out") (bold(v)))_2 <= 2 k (D - 2 r) (1 + ((D - 2 r)^2 / eta)^q)^(-m) -> 0 $
  as $D -> infinity$. Choosing $D_0$ large enough that $norm(nabla F_("out"))_2 < n_0 t_("min")^m r$ guarantees $chevron.l nabla F_(m, eta) (bold(v)), bold(u) chevron.r > 0$ everywhere on $partial B(bold(0), 2 r)$. By continuity of $F_(m, eta)$ on the compact ball $overline(B)(bold(0), 2 r)$, the minimum over $overline(B)(bold(0), 2 r)$ cannot occur on the boundary and is therefore an interior local minimizer.
]

#subsec("E", "High-Dimensional Scale Fragility and Mode Merging")

#proposition(title: "High-Dimensional Typicality Sensitivity", origin: "Original Derivation")[
  Suppose $bold(x)_i tilde.op cal(N)(bold(mu)_j, sigma^2 I_d)$ in dimension $d$ with inter-cluster separation $norm(bold(mu)_1 - bold(mu)_2)_2 = delta sqrt(d) sigma$, and let $eta_j = gamma d sigma^2$. By concentration of $chi_d^2$, within-cluster distances satisfy $d_(i j)^2 / eta_j -> 1 / gamma$ in probability as $d -> infinity$, while cross-cluster distances satisfy $d_(i l)^2 / eta_j -> (1 + delta^2) / gamma$. Consequently:
  1. If $gamma$ is underestimated ($gamma << 1$) with $m -> 1^+$, typicalities decay super-polynomially as $t_(i j) approx gamma^(1 / (m - 1)) -> 0$;
  2. If $delta^2$ is moderate relative to $gamma$, the two clusters' $eta$-spheres overlap heavily and the non-convex single-cluster functional $F_(m, eta)$ loses its secondary local basin, collapsing distinct prototypes into a single merged cluster.
]

#proof[
  For $bold(x) tilde.op cal(N)(bold(mu)_1, sigma^2 I_d)$, $norm(bold(x) - bold(mu)_1)_2^2 / (d sigma^2)$ has mean $1$ and variance $2 / d -> 0$, while $norm(bold(x) - bold(mu)_2)_2^2 / (d sigma^2)$ has mean $1 + delta^2$ and variance $(2 + 4 delta^2) / d -> 0$ @beyer1999nearest. Substituting these limits into Lemma 1 yields $t_(i 1) -> (1 + gamma^(-q))^(-1)$ and $t_(i 2) -> (1 + ((1 + delta^2) / gamma)^q)^(-1)$ in probability. When $gamma < 1$ and $q = 1 / (m - 1)$ is large, $t_(i 1) <= gamma^q -> 0$. Conversely, when $gamma$ is large enough that the $1 / 2$-typicality radius $sqrt(gamma d) sigma$ spans the inter-mode gap $delta sqrt(d) sigma$, the concave tail of $rho_(m, eta)$ blurs the inter-cluster valley, merging the two basins of $F_(m, eta)$ into one global mode centered at the mixture mean.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Because PCM updates each cluster independently, one AO sweep costs $O(n c d)$ arithmetic operations—identical asymptotically to FCM, and often faster in practice because scalar typicality evaluation $t_(i j) = (1 + (d_(i j)^2 / eta_j)^q)^(-1)$ avoids cross-cluster normalization loops and locks rapidly into local modes. Including the FCM warm-start ($T_0$ sweeps), total complexity is $O(n c d (T_0 + T))$ with memory $O(n c + n d + c d)$. Notably, once $bold(eta)$ is fixed, the $c$ clusters can run in completely independent parallel threads with zero synchronization.

Geometrically, PCM replaces FCM's Voronoi partition cells with $c$ overlapping radial typicality bells whose $1 / 2$-level sets are spheres (or $A$-ellipsoids) of radius $sqrt(eta_j)$ centered at $bold(v)_j$. Unlike FCM, PCM does not partition ambient space exhaustively: distant empty space has $max_j t_(i j) approx 0$, enabling simultaneous clustering and outlier detection.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n c d T)$ (+ $O(n c d T_0)$ FCM warm-start)],
    [Space Complexity], [$O(n c + n d + c d)$ (embarrassingly parallel in $c$)],
    [Inductive Geometry], [Radial $A$-spheres of scale $sqrt(eta_j)$ (Geman-McClure at $m = 2$)],
    [Hyperparameters], [Cluster count $c$, fuzzifier $m > 1$, bandwidth vector $bold(eta)$],
    [Robustness], [Strictly redescending influence; immune to distant noise],
    [Primary Pathology], [Coincident cluster collapse (Theorem 2)],
    [Convergence Type], [Monotone descent to local mode equilibria of $F_(m, eta_j)$],
  ),
  caption: [Operational and theoretical profile of Possibilistic C-Means (Alg. 152).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Noisy Medical Image Segmentation:* Ultrasound, low-dose CT, and MRI scans contain speckle artifacts and anatomical anomalies that violate closed-world tissue assumptions. PCM's redescending M-estimator weights (Theorem 1) suppress corrupted pixels automatically rather than letting them bias tissue intensity prototypes.
2. *Radar, Sonar, and LiDAR Target Extraction in Clutter:* Autonomous perception sensors collect dense background clutter alongside true objects. Because $t_(i j) -> 0$ smoothly outside each cluster's $eta_j$-neighborhood, clutter remains unassigned instead of distorting target centroids.
3. *Industrial Fault Monitoring and Open-Set Anomaly Detection:* In predictive maintenance, normal operating regimes form dense clusters while unprecedented faults appear far from known modes. The maximum typicality $max_j t_(i j)$ serves directly as a calibrated in-distribution compatibility score.
4. *Astronomical Source Separation from Diffuse Background:* Photometric catalogs mix compact stellar populations with diffuse sky background noise. PCM extracts dense stellar modes while rejecting unclustered background points without manual thresholding.
5. *Robust Document and Customer Prototyping with Spam:* Real web and transaction logs carry bot traffic and malformed records. Proposition 1 guarantees that even high fractions of distant spam records cannot pull local cluster prototypes away from genuine behavioral modes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

PCM's foremost failure mode is _coincident cluster collapse_ (Theorem 2), which strikes whenever clusters overlap substantially, differ strongly in density, or receive poorly calibrated $bold(eta)$. Second, _scale-parameter sensitivity_ causes under-estimated $eta_j$ to trap prototypes on isolated micro-clumps while over-estimated $eta_j$ merges distinct modes. Third, in _high-dimensional spaces_, Proposition 2 shows that distance concentration exacerbates mode merging, frequently reducing the effective unique cluster count to $1$ and rendering partition metrics undefined. Fourth, fixed Euclidean $A = I_d$ still misaligns with strongly elongated ellipsoids.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Cluttered Sensor Data], [Strong: redescending outlier weight], [Validate $bold(eta)$ via $alpha$-cuts],
    [Overlapping Clusters], [Poor: prototypes coalesce to one mode], [Use PFCM, FCM (Alg. 151), or GG],
    [High-$d$ Feature Spaces], [Prone to single-mode collapse], [Reduce $d$ first or use GK / GG],
    [Anisotropic Ellipsoids], [Poor under isotropic $A = I_d$], [Combine with GK covariance metric],
    [Open-Set Novelty], [Strong: $max_j t_(i j)$ flags outliers], [Compare against Student-$t$ mixtures @peel2000robust],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for PCM.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PCM sits at the intersection of fuzzy clustering, robust M-estimation, and density mode-seeking. Replacing FCM's simplex constraint with the penalty $eta_j (1 - t_(i j))^m$ converts a partitioner into $c$ decoupled Geman-McClure-type M-estimators (Theorem 1), closely paralleling trimmed $k$-means @cuesta1997trimmed, heavy-tailed Student-$t$ mixtures @peel2000robust, and kernel mode-seeking Mean Shift @cheng1995mean @comaniciu2002mean (whose fixed-point updates share the exact iteratively reweighted mean structure with decaying radial weights). To combine FCM's non-coincident partitioning with PCM's noise immunity, hybrid Possibilistic-Fuzzy C-Means (PFCM) optimizes coupled probabilistic and possibilistic terms jointly.

#proposition(title: "Monotone Typicality Defuzzification vs. Scale Bias", origin: "Original Proof")[
  Under PCM, assigning observation $bold(x)_i$ to $pi(i) in arg max_(1 <= j <= c) t_(i j)^*$ is equivalent to minimizing the _scale-normalized_ squared distance $d_(i j)^2 / eta_j$. Consequently, when $eta_1 = dots = eta_c$, hard defuzzification yields Euclidean Voronoi boundaries, whereas unequal $eta_j$ induces Apollonius hyperspherical boundaries favoring larger-scale clusters.
]

#proof[
  By Lemma 1, $t_(i j)^* = (1 + (d_(i j)^2 / eta_j)^(1 / (m - 1)))^(-1)$ is strictly decreasing in $z_(i j) = d_(i j)^2 / eta_j$. Hence $t_(i j)^* >= t_(i l)^*$ holds if and only if $d_(i j)^2 / eta_j <= d_(i l)^2 / eta_l$, whose equality locus $d_(i j) / d_(i l) = sqrt(eta_j / eta_l)$ is a hyperplane when $eta_j = eta_l$ and an Apollonius hypersphere enclosing the smaller-$eta$ prototype when $eta_j != eta_l$.
]

The ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing provides a textbook empirical confirmation of Theorems 1–2 and Proposition 2. On *Multi-Shape*, PCM achieves Silhouette $=$ *0.795* (ARI $=$ N/A)—the *highest* Silhouette among all four algorithms (Algs. 151–154), outperforming FCM ($0.713$), GK ($0.709$), and GG ($0.706$) because well-separated dense regions allow PCM's redescending weights to isolate clean modes while downweighting sparse transition points. Conversely, on overlapping datasets where basins compete, coincident drift degrades partition contrast: *Iris* drops to Silhouette $= 0.266$, ARI $= 0.510$ (versus FCM $0.447 / 0.594$), and *Customers* drops to Silhouette $= 0.188$ (versus FCM $0.316$). Most strikingly, on *High-D*, PCM reports Silhouette $=$ *N/A* (ARI $=$ N/A), directly manifesting Theorem 2 and Proposition 2: under high-dimensional distance concentration, decoupled modes merge into fewer than two distinct non-empty clusters, leaving Silhouette mathematically undefined. Finally, PCM's mean runtime of *0.627 ms* is the fastest of the quartet, reflecting its uncoupled scalar typicality updates and rapid local convergence.
