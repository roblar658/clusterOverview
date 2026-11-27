#import "../template/ieee.typ": *

#alg_entry(
  id: 156,
  name: "Spatial Fuzzy C-Means (sFCM)",
  category: "Cat. VII: Fuzzy, Biclustering & Domain",
  year: "2006",
  paradigm: "Spatially Regularized Fuzzy Partitioning",
  complexity: $O(T c n (d + |cal(N)|))$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Standard Fuzzy C-Means (FCM) @bezdek1981fcm treats observations
$bold(x)_1, dots, bold(x)_n in bb(R)^d$ as exchangeable and ignores their
spatial arrangement. In digital imaging, geospatial fields, and
sensor-array measurements, adjacent sites exhibit strong spatial
autocorrelation: neighboring pixels belong to the same homogeneous region
with high probability except along low-dimensional boundaries. When
additive Gaussian or impulsive salt-and-pepper noise corrupts individual
sites, pixelwise FCM produces speckled, topologically fragmented
partitions. Spatial Fuzzy C-Means (sFCM), introduced by Chuang et al.
@chuang2006sfcm, incorporates local neighborhood membership consensus
directly into the fuzzy membership update at every alternating iteration.

Let $G = (V, E)$ be an undirected spatial graph on $V = {1, dots, n}$
(for a 2D image lattice, $V$ indexes pixels and each site $i$ has a
closed local window $cal(N)_i subset.eq V$ of size $W_i = |cal(N)_i|$,
such as a $3 times 3$ or $5 times 5$ square neighborhood; for general
feature-space or sensor data, $cal(N)_i$ is the $k$-nearest-neighbor
spatial or metric neighborhood including $i$). Given cluster prototypes
$V_C = (bold(v)_1, dots, bold(v)_c) in (bb(R)^d)^c$ and fuzzifier
$m > 1$, the first stage of each iteration evaluates the classical
intensity-based FCM membership $u_(i j) in [0, 1]$ minimizing
$J_m (U, V_C) = sum_(i=1)^n sum_(j=1)^c u_(i j)^m norm(bold(x)_i - bold(v)_j)_2^2$:
$
  u_(i j) = ( sum_(k=1)^c
  (norm(bold(x)_i - bold(v)_j)_2^2 / norm(bold(x)_i - bold(v)_k)_2^2)^(1 / (m - 1)) )^(-1).
$
The second stage computes the _spatial consensus function_
$
  h_(i j) = sum_(l in cal(N)_i) u_(l j),
$
which aggregates the support for cluster $j$ across the neighborhood
$cal(N)_i$, and fuses feature and spatial evidence via parameterized
power weighting with exponents $p >= 0$, $q >= 0$ ($p + q > 0$):
$
  tilde(u)_(i j) = (u_(i j)^p h_(i j)^q) / (sum_(k=1)^c u_(i k)^p h_(i k)^q).
$
Cluster prototypes are then updated using the spatially regularized
memberships $tilde(u)_(i j)$:
$
  bold(v)_j^((t+1)) = (sum_(i=1)^n tilde(u)_(i j)^m bold(x)_i)
  / (sum_(i=1)^n tilde(u)_(i j)^m).
$

#definition(title: "Spatial Fuzzy C-Means Operator sFCM_(p,q)")[
  For a spatial neighborhood system ${cal(N)_i}_(i=1)^n$, fuzzifier
  $m > 1$, and exponents $(p, q) in bb(R)_(>= 0)^2 without {(0, 0)}$, the
  _sFCM spatial fusion map_ $cal(S)_(p,q) : cal(M)_(f c) -> cal(M)_(f c)$
  sends raw memberships $U$ to $tilde(U) = cal(S)_(p,q) (U)$ via
  $tilde(u)_(i j) = u_(i j)^p h_(i j)^q / sum_(k=1)^c u_(i k)^p h_(i k)^q$.
  Denoting the FCM membership update by $cal(F)_m (V_C)$ and the
  prototype update by $cal(G)_m (tilde(U))$, one iteration of
  $"sFCM"_(p,q)$ applies the composite fixed-point operator
  $cal(T)_(p,q) (V_C) = cal(G)_m (cal(S)_(p,q) (cal(F)_m (V_C)))$.
  Setting $(p, q) = (1, 0)$ recovers classical FCM @bezdek1981fcm.
]

#algbox(
  title: "Spatial Fuzzy C-Means (sFCM_(p,q))",
  inputs: [Data $X in bb(R)^(n times d)$, neighborhood sets ${cal(N)_i}_(i=1)^n$, cluster count $c$, fuzzifier $m > 1$, exponents $p, q >= 0$, tolerance $epsilon > 0$, max iterations $T$.],
  outputs: [Spatially smoothed membership matrix $tilde(U) in cal(M)_(f c)$ and prototypes $V_C in bb(R)^(c times d)$.],
)[
  1. Initialize prototypes $V_C^((0)) in bb(R)^(c times d)$ (or memberships $tilde(U)^((0)) in cal(M)_(f c)$).
  2. For $t = 0, 1, dots, T - 1$:
  3. *Feature-Space Membership Step:* For all $i, j$, compute $D_(i j)^2 = norm(bold(x)_i - bold(v)_j^((t)))_2^2$ and raw membership $u_(i j) = [sum_(k=1)^c (D_(i j)^2 / D_(i k)^2)^(1 / (m - 1))]^(-1)$ (splitting uniformly if $D_(i j)^2 = 0$).
  4. *Spatial Consensus Step:* For all $i = 1, dots, n$ and $j = 1, dots, c$, compute $h_(i j) = sum_(l in cal(N)_i) u_(l j)$.
  5. *Multiplicative Fusion Step:* Set $tilde(u)_(i j)^((t+1)) = u_(i j)^p h_(i j)^q / sum_(k=1)^c u_(i k)^p h_(i k)^q$ (preserving $u_(i j)$ if the denominator vanishes).
  6. *Prototype Update Step:* For each $j = 1, dots, c$, compute $bold(v)_j^((t+1)) = (sum_(i=1)^n (tilde(u)_(i j)^((t+1)))^m bold(x)_i) / (sum_(i=1)^n (tilde(u)_(i j)^((t+1)))^m)$.
  7. *Termination:* Stop if $norm(tilde(U)^((t+1)) - tilde(U)^((t)))_F < epsilon$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

Although sFCM is widely adopted for noisy image segmentation
@chuang2006sfcm, its multiplicative spatial operator introduces subtle
mathematical phenomena that differ fundamentally from classical
alternating minimization. Below we prove its exact mean-field Markov
Random Field (MRF) posterior structure, derive exact impulse-noise
rejection bounds, prove via an explicit counterexample that $"sFCM"_(p,q)$
is *not* a descent method for $J_m$ (contrasting with energy-penalized
spatial FCM), and quantify its critical thin-structure erasure threshold.

#subsec("A", "Invariance Properties and Mean-Field MRF Posterior Factorization")

#proposition(title: "Simplex Preservation, Scale Invariance, and MRF Posterior Form", origin: "Original Derivation")[
  Let $p, q >= 0$ with $p + q > 0$, and assume $u_(i j) > 0$ for all
  $i, j$. Then the spatial fusion operator $cal(S)_(p,q)$ satisfies:
  1. *(Simplex and Normalization Invariance)* $tilde(U) = cal(S)_(p,q) (U)$
     lies strictly in the open simplex product $(Delta_(c-1)^circle.small)^n$,
     and replacing $h_(i j)$ with the window average
     $bar(h)_(i j) = W_i^(-1) sum_(l in cal(N)_i) u_(l j)$ leaves
     $tilde(U)$ invariant.
  2. *(Affine Similarity Invariance)* Under uniform intensity scaling,
     rotation, and translation $bold(x)_i |-> lambda bold(Q) bold(x)_i + bold(b)$
     ($lambda > 0$, $bold(Q)^top bold(Q) = bold(I)_d$), the iterates
     $tilde(U)^((t))$ are exactly invariant and prototypes transform
     equivariantly as $bold(v)_j^((t)) |-> lambda bold(Q) bold(v)_j^((t)) + bold(b)$.
  3. *(Mean-Field MRF Factorization)* Defining likelihood
     $L_(i j) = norm(bold(x)_i - bold(v)_j)_2^(-2 / (m - 1))$ and
     spatial clique prior $pi_(i j) = bar(h)_(i j)^q$, the fused
     membership satisfies the exact tempered Bayes-MRF identity
     $
       tilde(u)_(i j) = (L_(i j)^p pi_(i j)) / (sum_(k=1)^c L_(i k)^p pi_(i k)).
     $
]

#proof[
  For (1), strict positivity of $u_(i j)$ implies $h_(i j) > 0$, hence
  $u_(i j)^p h_(i j)^q > 0$ and $sum_(j=1)^c tilde(u)_(i j) = 1$, placing
  each row in $Delta_(c-1)^circle.small$. Substituting $h_(i j) = W_i bar(h)_(i j)$
  factors $W_i^q > 0$ out of both numerator and denominator, canceling
  identically. For (2), under $bold(x)'_i = lambda bold(Q) bold(x)_i + bold(b)$
  and $bold(v)'_j = lambda bold(Q) bold(v)_j + bold(b)$, squared
  distances scale uniformly by $lambda^2$, which cancels in every ratio
  $D_(i j)^2 / D_(i k)^2$; hence $u_(i j)$, $h_(i j)$, and
  $tilde(u)_(i j)$ remain unchanged, and linearity of weighted averages
  propagates equivariance to $bold(v)'_j$. For (3), writing
  $u_(i j) = Z_i^(-1) L_(i j)$ with $Z_i = sum_(r=1)^c L_(i r)$ gives
  $u_(i j)^p bar(h)_(i j)^q = Z_i^(-p) L_(i j)^p pi_(i j)$. Canceling the
  cluster-independent factor $Z_i^(-p)$ between numerator and denominator
  yields the exact Bayes-MRF factorization @bishop2006pattern, showing
  that $p$ tempers the data log-likelihood while $q$ scales the
  logarithmic spatial consensus prior $log bar(h)_(i j)$.
]

#subsec("B", "Exact Impulse-Noise Rejection Threshold")

#theorem(title: "Impulse-Noise Restoration Guarantee in Homogeneous Regions", origin: "Original Proof")[
  Let pixel $i$ lie in a homogeneous region of true class $j^ast$ with
  neighborhood $cal(N)_i$ of size $W = |cal(N)_i| >= 2$ (containing $i$).
  Suppose the $W - 1$ uncorrupted neighbors satisfy $u_(l j^ast) >= 1 - eta$
  ($eta in [0, 1/2)$), while pixel $i$ is corrupted so that a rival
  class $k != j^ast$ dominates its raw membership with $u_(i k) = 1 - epsilon$
  and $u_(i j^ast) = epsilon in (0, 1/2)$. Then pixelwise FCM ($q = 0$)
  misclassifies pixel $i$, whereas $"sFCM"_(p,q)$ strictly restores
  $tilde(u)_(i j^ast) > tilde(u)_(i k)$ whenever
  $
    (epsilon / (1 - epsilon))^p
    (((W - 1)(1 - eta) + epsilon) / ((W - 1) eta + 1 - epsilon))^q > 1,
  $
  which is both necessary and sufficient when $c = 2$. In particular, for
  a clean background ($eta = 0$), $"sFCM"_(1,1)$ restores every impulse
  with $epsilon > 1 / (W + 1)$ (exact equivalence), and general
  $"sFCM"_(p,q)$ restores every impulse with
  $epsilon > (1 + (W - 1)^(q / p))^(-1)$.
]

#proof[
  For any rival cluster $k != j^ast$ and any uncorrupted neighbor
  $l in cal(N)_i without {i}$, simplex normalization gives
  $u_(l k) <= 1 - u_(l j^ast) <= eta$ (with equality when $c = 2$).
  Summing over $cal(N)_i$ therefore yields the tight bounds
  $
    h_(i j^ast) >= (W - 1)(1 - eta) + epsilon,
    quad h_(i k) <= (W - 1) eta + 1 - epsilon.
  $
  Since $tilde(u)_(i j^ast) > tilde(u)_(i k)$ iff
  $(u_(i j^ast) / u_(i k))^p (h_(i j^ast) / h_(i k))^q > 1$, substituting
  $u_(i j^ast) = epsilon$, $u_(i k) = 1 - epsilon$, and the bounds on
  $h_(i j^ast), h_(i k)$ yields the stated inequality (exact when $c = 2$).
  Setting $eta = 0$ and $p = q = 1$, the condition becomes
  $epsilon (W - 1 + epsilon) > (1 - epsilon)^2 = 1 - 2 epsilon + epsilon^2$.
  Canceling $epsilon^2$ gives $(W + 1) epsilon > 1$, i.e.
  $epsilon > 1 / (W + 1)$. For general $p, q > 0$ at $eta = 0$, observe
  that $(W - 1 + epsilon) / (1 - epsilon) > W - 1$; hence
  $epsilon / (1 - epsilon) > (W - 1)^(-q / p)$—equivalently
  $epsilon > (1 + (W - 1)^(q / p))^(-1)$—guarantees
  $(u_(i j^ast) / u_(i k))^p (h_(i j^ast) / h_(i k))^q > 1$.
]

For a standard $3 times 3$ window ($W = 9$), Theorem 1 proves that
$"sFCM"_(1,1)$ eliminates any isolated impulse outlier whose raw true-class
membership remains above $0.10$ ($9:1$ adverse intensity odds), while
$"sFCM"_(1,2)$ corrects impulses down to $epsilon > 1/65 approx 0.015$.

#subsec("C", "Non-Monotonicity of Classical J_m vs. Regularized Descent")

A critical distinction separates Chuang's post-hoc multiplicative fusion
@chuang2006sfcm from variational spatially regularized FCM @zhang2004kfcm:
for $q > 0$, $"sFCM"_(p,q)$ does *not* generally decrease the classical
FCM objective $J_m$.

#theorem(title: "Constructive Non-Monotonicity of sFCM_(1,q) on J_m and Variational Remedy", origin: "Original Counterexample")[
  1. *(Counterexample to $J_m$-Descent)* On the 1D lattice dataset
     $X = {0, 1, 2} subset.eq bb(R)$ with $c = 2$, $m = 2$, and 1-hop
     neighborhoods $cal(N)_1 = {1, 2}$, $cal(N)_2 = {1, 2, 3}$,
     $cal(N)_3 = {2, 3}$, initialized at the symmetric stationary
     minimizer $(U^ast, V_C^ast)$ of $J_2$, a single full iteration of
     $"sFCM"_(1,q)$ for *any* $q > 0$ strictly increases $J_2$:
     $J_2 (tilde(U), tilde(V)_C) > J_2 (U^ast, V_C^ast)$.
  2. *(Lyapunov-Monotone Variational Surrogate)* Replacing post-hoc
     fusion with the neighborhood-penalized energy @zhang2004kfcm
     $
       J_m^"reg" (U, V_C) = sum_(i=1)^n sum_(j=1)^c u_(i j)^m
       (norm(bold(x)_i - bold(v)_j)_2^2
       + alpha / W_i sum_(l in cal(N)_i) norm(bold(x)_l - bold(v)_j)_2^2)
     $
     ($alpha >= 0$) restores strict Lyapunov descent and Zangwill
     convergence @bezdek1981fcm @hathaway1986fcm.
]

#proof[
  For (1), parameterize symmetric prototypes on $X = {0, 1, 2}$ as
  $v_1 = 1 - a$, $v_2 = 1 + a$ ($a > 0$). For fixed $a$, exact $J_2$
  minimization over $U$ gives $u_(2 1) = u_(2 2) = 1/2$ and
  $u_(1 1) = u_(3 2) = (1 + a)^2 / (2(1 + a^2))$,
  $u_(1 2) = u_(3 1) = (1 - a)^2 / (2(1 + a^2))$. Using
  $min_(u in [0,1]) (u^2 d_1^2 + (1 - u)^2 d_2^2) = d_1^2 d_2^2 / (d_1^2 + d_2^2)$
  at each site yields the concentrated objective as a function of
  $s = a^2$:
  $
    g(s) = min_U J_2 (U, V(a)) = s / 2 + (1 - s)^2 / (1 + s)
    = 3 / 2 s - 3 + 4 / (1 + s).
  $
  Differentiating gives $g'(s) = 3/2 - 4 / (1 + s)^2$ and
  $g''(s) = 8 / (1 + s)^3 > 0$, so $g$ is strictly convex on $[0, infinity)$
  with unique minimizer $s^ast = sqrt(8/3) - 1 approx 0.6330$,
  $a^ast = sqrt(s^ast) approx 0.7956$. Checking the prototype update at
  $a^ast$ gives
  $1 - v_1^ast = (u_(1 1)^2 - u_(1 2)^2) / (u_(1 1)^2 + u_(1 2)^2 + 1/4) = 8 a (1 + a^2) / (3 + 14 a^2 + 3 a^4) = a^ast$
  since $3 (s^ast + 1)^2 = 8$; hence $(U^ast, V_C^ast)$ is an exact FCM
  stationary minimizer with $u_(1 1)^ast > 1/2 > u_(1 2)^ast$.

  Now execute one iteration of $"sFCM"_(1,q)$ ($q > 0$) from
  $(U^ast, V_C^ast)$. At site $1$,
  $h_(1 1) = u_(1 1)^ast + 1/2 > 1 > u_(1 2)^ast + 1/2 = h_(1 2)$, so
  $tilde(u)_(1 1) = u_(1 1)^ast h_(1 1)^q / (u_(1 1)^ast h_(1 1)^q + u_(1 2)^ast h_(1 2)^q) > u_(1 1)^ast$
  strictly. At site $2$, $h_(2 1) = u_(1 1)^ast + 1/2 + u_(1 2)^ast = 3/2 = h_(2 2)$,
  so $tilde(u)_(2 1) = 1/2$, and by symmetry $tilde(u)_(3 2) = tilde(u)_(1 1)$.
  The updated prototypes are $tilde(v)_1 = 1 - tilde(a)$,
  $tilde(v)_2 = 1 + tilde(a)$ with
  $tilde(a) = phi(tilde(u)_(1 1)) := (2 tilde(u)_(1 1) - 1) / (2 tilde(u)_(1 1)^2 - 2 tilde(u)_(1 1) + 5/4)$.
  Because $phi'(u) = (4 u (1 - u) + 1/2) / (2 u^2 - 2 u + 5/4)^2 > 0$
  on $[1/2, 1]$ and $phi(u_(1 1)^ast) = a^ast$, we obtain
  $tilde(a) > a^ast$ strictly. Since $g'(s) > 0$ for all $s > s^ast$,
  $
    J_2 (tilde(U), tilde(V)_C) >= min_U J_2 (U, tilde(V)_C)
    = g(tilde(a)^2) > g((a^ast)^2) = J_2 (U^ast, V_C^ast).
  $

  For (2), define effective distance
  $tilde(D)_(i j)^2 = norm(bold(x)_i - bold(v)_j)_2^2 + (alpha / W_i) sum_(l in cal(N)_i) norm(bold(x)_l - bold(v)_j)_2^2$.
  For fixed $V_C$, $U |-> J_m^"reg" (U, V_C)$ is strictly convex on
  $cal(M)_(f c)$ with unique closed-form minimizer using $tilde(D)_(i j)^2$.
  For fixed $U$, $V_C |-> J_m^"reg" (U, V_C)$ is strictly convex quadratic
  with unique minimizer
  $bold(v)_j = (sum_i u_(i j)^m (bold(x)_i + alpha bar(bold(x))_(cal(N)_i))) / ((1 + alpha) sum_i u_(i j)^m)$
  where $bar(bold(x))_(cal(N)_i) = W_i^(-1) sum_(l in cal(N)_i) bold(x)_l$.
  Alternating updates therefore decrease $J_m^"reg"$ monotonically and
  converge to stationary points by Zangwill's theorem @bezdek1981fcm
  @hathaway1986fcm.
]

#subsec("D", "Critical Thin-Structure Erasure Threshold")

Strong spatial weighting suppresses isolated noise at the cost of
eroding fine geometric structures. We derive the exact threshold at which
a 1-pixel linear anatomical feature is erased in a single spatial pass.

#proposition(title: "Exact Erasure Threshold for 1-Pixel Linear Structures", origin: "Original Counterexample")[
  Consider a 1-pixel-wide straight structure of class $2$ embedded inside
  a clean background of class $1$ ($u_(l 1) = 1$, $u_(l 2) = 0$) on a 2D
  grid, segmented with a $(2 r + 1) times (2 r + 1)$ square window
  ($r >= 1$). Suppose every pixel $i$ on the thin structure has raw
  feature membership $u_(i 2) = 1 - epsilon$ and $u_(i 1) = epsilon in (0, 1/2)$.
  Then a single spatial fusion pass of $"sFCM"_(p,q)$ erases the
  structure ($tilde(u)_(i 1) > tilde(u)_(i 2)$) if and only if
  $
    ((1 - epsilon) / epsilon)^p < ((2 r + epsilon) / (1 - epsilon))^q.
  $
  In particular, for a $5 times 5$ window ($r = 2$) under $"sFCM"_(1,2)$
  @chuang2006sfcm, even a high-confidence feature with $u_(i 2) = 0.94$
  ($epsilon = 0.06$) is completely erased in one iteration.
]

#proof[
  Inside a $(2 r + 1) times (2 r + 1)$ window centered on a pixel $i$ of
  the straight 1-pixel line, exactly $2 r + 1$ sites lie on the line and
  $(2 r + 1)^2 - (2 r + 1) = 2 r (2 r + 1)$ sites belong to the clean
  background. Summing raw memberships over $cal(N)_i$ gives
  $
    h_(i 2) = (2 r + 1)(1 - epsilon), \
    h_(i 1) = 2 r (2 r + 1) + (2 r + 1) epsilon = (2 r + 1)(2 r + epsilon).
  $
  The common factor $2 r + 1$ cancels identically in the spatial ratio
  $h_(i 1) / h_(i 2) = (2 r + epsilon) / (1 - epsilon)$. Erasure occurs
  ($tilde(u)_(i 1) > tilde(u)_(i 2)$) iff
  $u_(i 1)^p h_(i 1)^q > u_(i 2)^p h_(i 2)^q$, which rearranges directly
  to $((1 - epsilon) / epsilon)^p < ((2 r + epsilon) / (1 - epsilon))^q$.
  Evaluating at $r = 2$, $p = 1$, $q = 2$, and $epsilon = 0.06$ yields
  LHS $= 0.94 / 0.06 approx 15.67 < (4.06 / 0.94)^2 approx 18.65 =$ RHS,
  confirming immediate single-step erasure.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of $"sFCM"_(p,q)$ executes four linear passes over the
$n$ sites: (i) squared distances and raw FCM memberships in $O(c n d)$
time; (ii) spatial neighborhood aggregation $h_(i j) = sum_(l in cal(N)_i) u_(l j)$
in $O(c n |cal(N)|)$ time; (iii) multiplicative fusion in $O(c n)$ time;
and (iv) prototype updates in $O(c n d)$ time. Across $T$ iterations,
total runtime is $O(T c n (d + |cal(N)|))$ with $O(n (d + c))$ memory
(plus $O(n |cal(N)|)$ if explicit neighbor lists are stored for
non-grid graphs). On regular 2D/3D image lattices, the spatial window sum
$h_(i j)$ is a separable box convolution that reduces from
$O(c n (2 r + 1)^2)$ to $O(c n)$ independent of window radius $r$ via
integral images (summed-area tables) or 1D sliding sums! For general
tabular or point-cloud datasets without an intrinsic grid, precomputing
a $k$-NN spatial graph adds a one-time $O(n^2 d)$ or $O(n d log n)$
setup cost.

Geometrically, sFCM lifts pixelwise Voronoi cells into a coupled
product of intensity space and spatial topology: homogeneous interior
regions experience strong membership sharpening (entropy reduction),
while straight interfaces remain stable when $p >= q$ and high-curvature
corners or thin filaments shrink inward according to Proposition 2.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(T c n (d + |cal(N)|))$; $O(T c n d)$ on grids via integral images],
    [Space Complexity], [$O(n (d + c))$ on regular grids; $O(n (d + c + |cal(N)|))$ on graphs],
    [Spatial Coupling], [Multiplicative log-consensus MRF prior $tilde(u)_(i j) prop u_(i j)^p h_(i j)^q$],
    [Noise Tolerance], [Exact impulse rejection for $epsilon > (1 + (W - 1)^(q / p))^(-1)$],
    [Hyperparameters], [Cluster count $c$, fuzzifier $m > 1$, window $cal(N)_i$, exponents $(p, q)$],
    [Convergence], [Fixed-point iteration (not $J_m$-monotone for $q > 0$); fast empirical stability],
  ),
  caption: [Operational, memory, and theoretical profile of Spatial Fuzzy C-Means (sFCM).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Neuroimaging MRI Tissue Classification:* T1/T2/PD brain scans
   contain spatially contiguous anatomical tissues (CSF, gray matter,
   white matter) corrupted by thermal noise and RF inhomogeneity.
   Chuang et al. @chuang2006sfcm demonstrated that $"sFCM"_(1,1)$ and
   $"sFCM"_(1,2)$ dramatically reduce misclassification salt-and-pepper
   artifacts while preserving soft partial-volume transitions across
   tissue interfaces.
2. *Ultrasound and Echocardiography Speckle Segmentation:* B-mode
   ultrasound exhibits multiplicative Rayleigh speckle whose pixelwise
   intensity distributions heavily overlap between cardiac chambers and
   myocardium. Aggregating memberships over local spatial windows
   suppresses speckle variance by a factor scaling with $W$ (Theorem 1),
   enabling coherent endocardial boundary extraction.
3. *Multispectral Satellite Land-Cover Mapping:* High-resolution
   earth-observation rasters suffer single-pixel cloud haze, sensor
   striping, and mixed-pixel canopy variability. Applying sFCM on 2D
   pixel lattices enforces spatial parcel coherence across agricultural
   fields and water bodies while retaining fuzzy boundary grades for
   ecotone transitions.
4. *Digital Pathology and Histological Stain Quantification:* Whole-slide
   biopsy tiles exhibit subcellular texture noise within macroscopically
   contiguous epithelial, stromal, and necrotic regions. Moderate spatial
   weighting ($p = 1, q = 1$) averages over local staining heterogeneity
   to segment coherent tissue compartments without discarding soft
   infiltration zones.
5. *Geophysical Subsurface and Spatial Transcriptomics Domains:* Spatial
   transcriptomics spots and seismic tomography grids associate each
   high-dimensional feature vector with physical 2D/3D coordinates.
   Defining $cal(N)_i$ via physical adjacency graphs allows sFCM to
   recover spatially coherent anatomical niches or lithological strata
   even when individual spot read counts are sparse and noisy.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four failure regimes govern practical deployment. First, *thin-structure
and small-lesion erasure* (Proposition 2) occurs whenever window radius
$r$ or spatial exponent $q / p$ exceeds the contrast-to-width ratio of
blood vessels, cortical sulci, or micro-calcifications. Second, *absence
of spatial homophily* causes sFCM to degrade accuracy when neighboring
sites belong to distinct classes (high-frequency checkerboard textures or
unstructured tabular records where $k$-NN neighborhoods span class
boundaries). Third, *non-monotone boundary oscillation* (Theorem 2) can
prevent strict $J_m$ convergence under aggressive $q >> p$, occasionally
causing small limit cycles on sharp edges unless damping or early
stopping is applied. Fourth, *non-convex feature clusters* remain
inseparable because the underlying feature distance is still Euclidean.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Noisy Medical / Sonar Images], [Excellent: strong impulse & Gaussian denoising], [Keep $3 times 3$ window and $q <= p$ near thin vessels],
    [Fine Angiography / Thin Filaments], [Poor for large $r, q$: erases 1-pixel structures], [Use anisotropic diffusion, Watershed @vincent1991watersheds, or small $q$],
    [Superpixel / Edge-Aligned Vision], [Moderate: isotropic windows blur sharp edges], [Use SLIC @achanta2012slic or Normalized Cut @shi2000ncut],
    [Curved Feature Manifolds + Noise], [Moderate: Euclidean feature step is linear], [Combine with RKHS kernel via KFCM-S @zhang2004kfcm (Alg. 155)],
    [Unstructured Tabular Data], [Context-dependent: relies on $k$-NN homophily], [Verify neighborhood purity; compare standard FCM @bezdek1981fcm],
  ),
  caption: [Domain suitability and failure-mode matrix for Spatial Fuzzy C-Means (sFCM).],
)

#sec("V", "Structural Relations and Empirical Behavior")

Spatial FCM occupies a pivotal bridge between classical prototype-based
fuzzy clustering and graph-regularized spatial segmentation. Formally, we
characterize its exact nested relationships across parameter limits.

#proposition(title: "Asymptotic Reductions of sFCM_(p,q) Across Exponent Regimes", origin: "Original Derivation")[
  Consider the operator $cal(T)_(p,q)$ of $"sFCM"_(p,q)$ on a connected
  neighborhood graph $G = (V, E)$ with $u_(i j) > 0$ everywhere:
  1. *(FCM Limit)* Setting $(p, q) = (1, 0)$ yields exact classical FCM
     @bezdek1981fcm; more generally, $(p, 0)$ with $p > 0$ is algebraically
     identical to classical FCM with effective fuzzifier
     $m_"eff" = 1 + (m - 1) / p$.
  2. *(Pure Consensus Diffusion Limit)* Setting $(p, q) = (0, 1)$ reduces
     $cal(S)_(0,1)$ to row-normalized graph random-walk averaging
     $tilde(U) = bold(D)^(-1) bold(A)_G U$, where $bold(A)_G$ is the
     adjacency matrix with self-loops and $bold(D) = "diag"(W_1, dots, W_n)$.
  3. *(Hard Local Majority Limit)* As $q -> infinity$ with $p$ fixed,
     $tilde(u)_(i j) -> bb(I)(j = arg max_k h_(i k))$ whenever the local
     neighborhood plurality is unique, implementing a crisp Potts/MRF
     majority filter.
]

#proof[
  For (1), when $q = 0$, $tilde(u)_(i j) = u_(i j)^p / sum_k u_(i k)^p$.
  Substituting $u_(i j) prop (D_(i j)^2)^(-1 / (m - 1))$ gives
  $tilde(u)_(i j) prop (D_(i j)^2)^(-p / (m - 1)) = (D_(i j)^2)^(-1 / (m_"eff" - 1))$,
  where $m_"eff" - 1 = (m - 1) / p > 0$, matching standard FCM with
  exponent $m_"eff"$ identically. For (2), when $p = 0, q = 1$,
  $sum_(k=1)^c h_(i k) = sum_(l in cal(N)_i) sum_(k=1)^c u_(l k) = W_i$,
  so $tilde(u)_(i j) = W_i^(-1) sum_(l in cal(N)_i) u_(l j) = (bold(D)^(-1) bold(A)_G U)_(i j)$,
  a single step of lazy random-walk diffusion on $G$. For (3), let
  $j^ast = arg max_k h_(i k)$ be unique so $h_(i k) / h_(i j^ast) < 1$
  for all $k != j^ast$. Dividing numerator and denominator of
  $tilde(u)_(i j)$ by $h_(i j^ast)^q$ drives all terms $k != j^ast$ to
  $0$ exponentially as $q -> infinity$, leaving $tilde(u)_(i j^ast) -> 1$.
]

When $k$-NN feature graphs supply $cal(N)_i$ on non-image data, sFCM acts
similarly to manifold diffusion and label-propagation regularization,
contrasting with KFCM (Alg. 155) @zhang2004kfcm which alters the global
metric via a Mercer kernel, and with Possibilistic C-Means
@krishnapuram1993pcm which achieves robustness by relaxing simplex
normalization instead of leveraging spatial context.

*Empirical Benchmark Behavior (ClusterOverview Suite).* Across the
four-dataset ClusterOverview benchmark (evaluated via Silhouette
@rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing),
sFCM achieves the strongest overall geometric cohesion among Algorithms
155–158. On *Multi-Shape*, it ties KFCM at $"Sil" = 0.713$
($"ARI" = "N/A"$), showing that neighborhood membership diffusion
propagates cluster coherence cleanly along locally smooth manifolds. On
*Customers* ($"Sil" = 0.370$) and *High-D* ($"Sil" = 0.303$,
$"ARI" = "N/A"$), sFCM substantially outperforms KFCM ($0.316$ and
$0.185$), Cheng-Church ($0.116$ and $0.043$), and PLAID ($0.039$ and
$0.072$): local neighborhood consensus stabilizes noisy assignments in
high-dimensional space without suffering from global Gaussian kernel
bandwidth collapse. On *Iris*, sFCM attains $"Sil" = 0.449$ (marginally
above KFCM's $0.447$) yet slightly lower ground-truth alignment
($"ARI" = 0.563$ vs. KFCM's $0.594$), directly reflecting the boundary
smoothing effect of Proposition 2 across the closely abutting
Versicolor–Virginica transition. Its mean runtime of $2.409 "ms"$ is
$11.5\%$ faster than KFCM ($2.723 "ms"$), consistent with its linear
$O(T c n (d + |cal(N)|))$ scaling once local neighborhoods are indexed.
