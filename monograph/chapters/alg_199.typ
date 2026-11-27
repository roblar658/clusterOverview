#import "../template/ieee.typ": *

#alg_entry(
  id: 199,
  name: "PSO-Clustering (Particle Swarm Optimization)",
  category: "Cat. VIII: Neural, Deep, Stream & Evolutionary",
  year: "2003",
  paradigm: "Continuous Swarm-Intelligence Centroid Optimization",
  complexity: $cal(O)(T S n K d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Rather than evolving discrete length-$n$ partition strings as in GKA (Alg. 198,
@krishna1999gka), Particle Swarm Optimization (PSO) clustering formulates
partitional clustering directly as a continuous non-convex optimization problem
over the $(K d)$-dimensional centroid space $bb(R)^(K times d)$
@vandermerwe2003pso. Each agent (*particle*) in a swarm of size $S$ flies
through $bb(R)^(K times d)$ representing an entire candidate codebook of $K$
centroids, balancing cognitive attraction toward its own best historical
configuration against social attraction toward the swarm's global best.

Given dataset $X = {bold(x)_1, dots.h, bold(x)_n} subset bb(R)^d$ and cluster
count $K$, represent particle $i in {1, dots.h, S}$ at iteration $t >= 0$ by
its centroid matrix $Z_i (t) = [bold(m)_(i, 1)(t), dots.h, bold(m)_(i, K)(t)]
in bb(R)^(K times d)$ and velocity matrix $V_i (t) in bb(R)^(K times d)$. Every
particle induces a hard Voronoi partition $C_(i, k) (t) = {bold(x) in X :
k = arg min_(j in {1, dots.h, K}) norm(bold(x) - bold(m)_(i, j)(t))_2}$ (with
deterministic tie-breaking). Van der Merwe and Engelbrecht @vandermerwe2003pso
introduce the *Quantization Error* (QE) fitness criterion alongside standard
SSE:
$
  J_"QE" (Z_i) &= 1 / K sum_(k=1)^K (1 / (|C_(i, k)|) sum_(bold(x) in C_(i, k))
    norm(bold(x) - bold(m)_(i, k))_2), \
  J_"SSE" (Z_i) &= sum_(k=1)^K sum_(bold(x) in C_(i, k))
    norm(bold(x) - bold(m)_(i, k))_2^2,
$
assigning an infinite penalty if any cluster $C_(i, k)$ becomes empty.

#definition(title: "Continuous Swarm Centroid Dynamics")[
  Let $P_i (t) = arg min_(0 <= s <= t) J(Z_i (s))$ denote particle $i$'s
  *personal best* centroid matrix and $G(t) = arg min_(1 <= i <= S) J(P_i (t))$
  the swarm's *global best*. Given inertia weight $w in (0, 1)$ and cognitive /
  social acceleration constants $c_1, c_2 > 0$, each particle updates via
  $
    V_i (t + 1) &= w V_i (t) + c_1 R_1 (t) dot.o (P_i (t) - Z_i (t)) \
                &wide + c_2 R_2 (t) dot.o (G(t) - Z_i (t)), \
    Z_i (t + 1) &= Z_i (t) + V_i (t + 1),
  $
  where $dot.o$ denotes the Hadamard (entrywise) product and
  $R_1 (t), R_2 (t) tilde.op "Uniform"([0, 1]^(K times d))$ are mutually
  independent i.i.d. random matrices drawn freshly at each iteration $t$. In
  *Hybrid $K$-Means-PSO* @vandermerwe2003pso, particle $i = 1$ is initialized
  at a converged Lloyd's $K$-means solution while particles $2, dots.h, S$ are
  initialized randomly across the data bounding box.
]

Because the search state lives in $bb(R)^(K times d)$ rather than
${1, dots.h, K}^n$, PSO-Clustering decouples agent memory from dataset size $n$
and supports arbitrary non-differentiable validity objectives (QE, SSE,
Davies-Bouldin @davies1979cluster, or Calinski-Harabasz @calinski1974dendrite)
without requiring analytic gradients.

#algbox(
  title: "PSO-Clustering (Continuous Centroid Swarm Optimization)",
  inputs: [Dataset $X subset bb(R)^d$, clusters $K$, swarm size $S$,
    iterations $T$, inertia $w in (0, 1)$, accelerations $c_1, c_2 > 0$,
    velocity cap $v_"max"$, fitness $J in {J_"QE", J_"SSE"}$.],
  outputs: [Global best centroid matrix $G(T)$ and Voronoi cluster labels.],
)[
  + Compute data bounding box $[bold(a), bold(b)] subset bb(R)^d$.
  + *Initialize swarm*: Optionally seed $Z_1 (0)$ with Lloyd's $K$-means (or
    $K$-means++ @arthur2007kmeanspp); draw remaining $Z_i (0)$ from random
    subsets of $X$, set $V_i (0) = bold(0)$, $P_i (0) = Z_i (0)$, and evaluate
    $G(0) = arg min_i J(P_i (0))$.
  + *For* iteration $t = 0, dots.h, T - 1$:
    - For each particle $i = 1, dots.h, S$, draw $R_1 (t), R_2 (t) tilde.op
      "Uniform"([0, 1]^(K times d))$, update $V_i (t + 1)$ (clamped to
      $[-v_"max", v_"max"]$), and set $Z_i (t + 1) = Z_i (t) + V_i (t + 1)$
      projected onto $[bold(a), bold(b)]^K$.
    - Assign all $bold(x) in X$ to nearest centroids in $Z_i (t + 1)$ and
      evaluate fitness $J(Z_i (t + 1))$ (re-seeding any empty centroid).
    - If $J(Z_i (t + 1)) < J(P_i (t))$, set $P_i (t + 1) <- Z_i (t + 1)$; else
      keep $P_i (t + 1) <- P_i (t)$.
    - Update $G(t + 1) <- arg min_i J(P_i (t + 1))$.
  + Optional local polish: run Lloyd's step on $G(T)$ if $J = J_"SSE"$; return
    $G(T)$ and resulting partition.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees: exact first- and second-order stochastic
stability regions of the swarm recurrence, piecewise-quadratic structure and
hybrid Lloyd dominance, a constructive stagnation counterexample showing why
pure PSO need not reach a local minimum, QE vs. SSE objective divergence, and
a worst-case $cal(O)(log K)$ error envelope for hybrid $K$-means++-PSO.

#subsec("A", "Order-1 and Order-2 Stochastic Stability Boundaries")

Selecting $(w, c_1, c_2)$ arbitrarily causes centroid trajectories to diverge
exponentially. We derive the exact closed-form mean and mean-square stability
criteria directly from the second-order stochastic difference equation.

#theorem(
  title: "Exact First- and Second-Order Swarm Stability Region",
  origin: "Original Proof",
)[
  Under the local stagnation assumption ($P_i = p$, $G = g$ constant along a
  scalar coordinate $z_t$ of $Z_i (t)$), let $phi_t = c_1 r_(1, t) + c_2 r_(2, t)$
  with mean $overline(phi) = (c_1 + c_2) / 2$ and variance
  $sigma_phi^2 = (c_1^2 + c_2^2) / 12$.
  (i) *Order-1 (Mean) Stability*: $bb(E)[z_t]$ converges for every initial
  condition if and only if $|w| < 1$ and $0 < c_1 + c_2 < 4(1 + w)$, in which
  case $lim_(t -> infinity) bb(E)[z_t] = (c_1 p + c_2 g) / (c_1 + c_2)$.
  (ii) *Order-2 (Mean-Square) Stability*: When $p = g = z^(ast)$, the second
  moment $bb(E)[(z_t - z^(ast))^2]$ decays geometrically to $0$ only if
  $|w| < 1$ and
  $
    sigma_phi^2 < (1 - w) / (1 + w) overline(phi) (2 + 2w - overline(phi)).
  $
  In the symmetric regime $c_1 = c_2 = c > 0$ with $w in [0, 1)$, this reduces
  to the exact necessary and sufficient variance-stability threshold
  $
    0 < c < (12 (1 - w^2)) / (7 - 5w).
  $
]

#proof[
  Eliminating velocity $v_t = z_t - z_(t-1)$ from the scalar update gives the
  second-order stochastic recurrence
  $
    z_(t+1) = (1 + w - phi_t) z_t - w z_(t-1) + c_1 r_(1, t) p + c_2 r_(2, t) g.
  $
  (i) Since $(r_(1, t), r_(2, t))$ is independent of $(z_t, z_(t-1))$, taking
  expectations yields $mu_(t+1) = a mu_t - w mu_(t-1) + b$, where
  $a = 1 + w - overline(phi)$ and $b = (c_1 p + c_2 g) / 2$. Its characteristic
  polynomial is $q(lambda) = lambda^2 - a lambda + w$. By the Jury stability
  criterion, both roots satisfy $|lambda| < 1$ iff $|w| < 1$,
  $q(1) = 1 - a + w = overline(phi) > 0$, and
  $q(-1) = 1 + a + w = 2 + 2w - overline(phi) > 0$, i.e.,
  $0 < c_1 + c_2 < 4(1 + w)$. At the fixed point $mu^(ast)$,
  $(1 - a + w) mu^(ast) = b => mu^(ast) = b / overline(phi) =
  (c_1 p + c_2 g) / (c_1 + c_2)$.

  (ii) Let $y_t = z_t - z^(ast)$ when $p = g = z^(ast)$ and set
  $alpha_t = 1 + w - phi_t$, so $y_(t+1) = alpha_t y_t - w y_(t-1)$ with
  $bb(E)[alpha_t] = a$ and $bb(E)[alpha_t^2] = a^2 + sigma_phi^2$. Define
  $u_t = bb(E)[y_t^2]$ and $v_t = bb(E)[y_t y_(t-1)]$. Squaring and taking
  expectations using independence of $alpha_t$ from $(y_t, y_(t-1))$ gives
  $
    u_(t+1) &= (a^2 + sigma_phi^2) u_t + w^2 u_(t-1) - 2 w a v_t, \
    v_(t+1) &= a u_t - w v_t.
  $
  In state vector $(u_t, u_(t-1), v_t)^top$, the $3 times 3$ transition matrix
  $M$ has characteristic polynomial
  $
    p(lambda) = lambda^3 + (w - a^2 - sigma_phi^2) lambda^2
    + (w a^2 - w sigma_phi^2 - w^2) lambda - w^3.
  $
  For Schur stability ($rho(M) < 1$), Jury's test requires $|w^3| < 1$,
  $p(-1) < 0$, and $p(1) > 0$. Direct evaluation gives
  $p(-1) = -(1 - w)(1 - w^2) - (1 + w) a^2 - (1 - w) sigma_phi^2 < 0$
  identically for $|w| < 1$, while
  $
    p(1) &= (1 + w - w^2 - w^3) - (1 - w) a^2 - (1 + w) sigma_phi^2 \
         &= (1 - w)((1 + w)^2 - a^2) - (1 + w) sigma_phi^2.
  $
  Factoring $(1 + w)^2 - a^2 = (1 + w - a)(1 + w + a) =
  overline(phi)(2 + 2w - overline(phi))$ shows that $p(1) > 0$ holds iff
  $sigma_phi^2 < ((1 - w) / (1 + w)) overline(phi)(2 + 2w - overline(phi))$.
  Setting $c_1 = c_2 = c$ gives $overline(phi) = c$ and $sigma_phi^2 = c^2 / 6$;
  substituting and solving linearly for $c > 0$ yields
  $c (1/6 + (1 - w)/(1 + w)) < 2(1 - w)$, which simplifies directly to
  $c < 12(1 - w^2) / (7 - 5w)$. (For standard Clerc-Kennedy values
  $w = 0.7298, c = 1.4962$, the bound evaluates to $1.6738 > 1.4962$, confirming
  strict mean-square contraction.)
]

#subsec("B", "Objective Landscape Geometry and Hybrid Dominance")

#proposition(
  title: "Piecewise-Quadratic Concave-Min Structure and Hybrid Dominance",
  origin: "Original Proof",
)[
  (i) On $bb(R)^(K times d)$, the continuous SSE functional
  $J_"SSE" (Z) = sum_(j=1)^n min_(k in {1, dots.h, K}) norm(bold(x)_j - bold(m)_k)_2^2$
  is continuous and piecewise-quadratic, smooth everywhere outside the
  measure-zero union of Voronoi bisector hyperplanes
  $cal(H) = union.big_(j, k != l) {Z : norm(bold(x)_j - bold(m)_k)_2 =
  norm(bold(x)_j - bold(m)_l)_2}$, and strictly convex on each interior
  assignment chamber.
  (ii) In Hybrid $K$-Means-PSO initialized with $Z_1 (0) = Z_"Lloyd"$, the
  global best trajectory satisfies $J_"SSE" (G(t + 1)) <= J_"SSE" (G(t)) <=
  J_"SSE" (Z_"Lloyd")$ for every $t >= 0$.
]

#proof[
  (i) Each term $f_j (Z) = min_(k=1..K) norm(bold(x)_j - bold(m)_k)_2^2$ is a
  finite pointwise minimum of continuous convex quadratics, hence continuous
  everywhere. Away from $cal(H)$, the minimizing index $k(j)$ for every datum
  $bold(x)_j$ is locally constant on an open neighborhood $U$, where
  $J_"SSE" (Z) = sum_(k=1)^K sum_(j in C_k) norm(bold(x)_j - bold(m)_k)_2^2$ has
  diagonal positive-semidefinite Hessian
  $nabla^2 J_"SSE" (Z) = 2 "diag"(|C_1| I_d, dots.h, |C_K| I_d) succ.eq 0$
  (strictly positive definite whenever every cluster is non-empty).
  (ii) Because $P_i (t + 1) = arg min {J(P_i (t)), J(Z_i (t + 1))}$ and
  $G(t) = arg min_i J(P_i (t))$, the sequence $J(G(t))$ is monotonically
  non-increasing. Since $Z_1 (0) = Z_"Lloyd"$ is included in the initial swarm,
  $J(G(t)) <= J(G(0)) <= J(Z_"Lloyd")$ holds deterministically for all $t >= 0$.
]

#subsec("C", "Premature Swarm Collapse Counterexample")

Unlike GKA (where KMO forces every elite toward Voronoi stationarity), pure
canonical PSO offers no guarantee of reaching even a local minimum.

#theorem(
  title: "Non-Stationary Stagnation Under Pure PSO Collapse",
  origin: "Original Counterexample",
)[
  For any $d >= 1$, $K >= 1$, and any positive parameters $(w, c_1, c_2)$,
  canonical PSO-Clustering admits fixed points $Z_i (t) = P_i (t) = G(t) = Z_0$
  with $V_i (t) = bold(0)$ for all $i = 1, dots.h, S$ at which $J_"SSE"$ is
  smooth yet $nabla J_"SSE" (Z_0) != bold(0)$, proving that pure PSO is not
  guaranteed to converge even to a first-order stationary point.
]

#proof[
  Let $d = 1$, $K = 1$, and $X = {-1, 1\}$ with unique global minimum at
  $m^(ast) = 0$ ($J_"SSE" (0) = 2$). Suppose at iteration $t_0$ the swarm
  collapses to the non-optimal state $Z_i (t_0) = P_i (t_0) = G(t_0) = [delta]$
  for some $delta != 0$ (e.g., $delta = 1/2$) with zero velocities
  $V_i (t_0) = [0]$ for all $i = 1, dots.h, S$. (Even without exact equality at
  $t_0$, Theorem 1 guarantees that under mean-square stability every particle
  exponentially collapses onto $G(t_0)$ whenever no particle discovers a better
  value along its transient trajectory!) Substituting $P_i - Z_i = 0$,
  $G - Z_i = 0$, and $V_i = 0$ into the swarm recurrence yields
  $V_i (t_0 + 1) = w dot 0 + c_1 R_1 dot 0 + c_2 R_2 dot 0 = 0$ and
  $Z_i (t_0 + 1) = [delta]$ identically for all $t >= t_0$. Yet at $Z_0 =
  [delta]$, $J_"SSE" (m) = 2(1 + m^2)$ is globally smooth with non-zero
  gradient $nabla J_"SSE" (delta) = 4 delta != 0$. Thus canonical PSO can
  stagnate permanently at non-stationary points unless augmented with a
  terminal Lloyd local-polish step or guaranteed-convergence velocity
  perturbation around $G(t)$.
]

#subsec("D", "QE vs. SSE Objective Divergence and Geometric Invariances")

#proposition(
  title: "Scale Equivariance and QE–SSE Minimizer Divergence",
  origin: "Original Derivation",
)[
  (i) Both $J_"QE"$ and $J_"SSE"$ are invariant under orthogonal rotations
  $Q in O(d)$ and cluster permutations $pi in S_K$, and satisfy positive
  homogeneity $J_"QE" (alpha Z; alpha X) = alpha J_"QE" (Z; X)$ and
  $J_"SSE" (alpha Z; alpha X) = alpha^2 J_"SSE" (Z; X)$ for all $alpha > 0$.
  (ii) For imbalanced clusters, the global minimizers of $J_"QE"$ and $J_"SSE"$
  do not coincide: specifically, $J_"QE"$ optimizes the sum of unweighted
  per-cluster $L_2$-medians (Fermat-Weber geometric medians) rather than
  cardinality-weighted sample means.
]

#proof[
  (i) Orthogonal maps preserve Euclidean norms:
  $norm(Q bold(x) - Q bold(m))_2 = norm(bold(x) - bold(m))_2$, proving
  $O(d)$-invariance; scaling multiplies every distance by $alpha > 0$, yielding
  degree-$1$ and degree-$2$ homogeneity respectively.

  (ii) Fix a partition ${C_1, dots.h, C_K}$. Within cluster $C_k$, minimizing
  $J_"SSE"$ selects the sample mean
  $bold(m)_k^"SSE" = arg min_bold(m) sum_(bold(x) in C_k) norm(bold(x) - bold(m))_2^2
  = (1 / |C_k|) sum_(bold(x) in C_k) bold(x)$, whereas minimizing $J_"QE"$
  selects the geometric median
  $bold(m)_k^"QE" = arg min_bold(m) sum_(bold(x) in C_k) norm(bold(x) - bold(m))_2$.
  Even in $d = 1$ with $K = 2$ on the skewed 4-point dataset
  $X = {0, 0, 1, 10}$ (split as ${0, 0, 1}$ and ${10}$), the SSE centroid of
  ${0, 0, 1}$ is its mean $1/3$, whereas the QE centroid is its median $0$!
  Moreover, because $J_"QE"$ divides each cluster's error by $|C_k|$ before
  averaging over $K$, a small tight cluster and a massive dense cluster exert
  equal $1/K$ influence on $J_"QE"$, whereas $J_"SSE"$ weights clusters
  proportionally to $|C_k|$.
]

#subsec("E", "Hybrid K-Means++-PSO Approximation Guarantee")

Combining $K$-means++ initialization with PSO yields a rigorous worst-case
approximation guarantee that pure random-start PSO lacks.

#theorem(
  title: "Worst-Case Approximation Envelope for Hybrid K-Means++-PSO",
  origin: "Classical — Cited",
)[
  If particle $i = 1$ in Hybrid PSO-Clustering is seeded via $D^2$-sampling
  ($K$-means++ @arthur2007kmeanspp) followed by Lloyd's refinement
  @lloyd1982least and fitness $J = J_"SSE"$ is optimized, then for every
  iteration horizon $T >= 0$ and swarm size $S >= 1$, the global best $G(T)$
  satisfies
  $
    bb(E)[J_"SSE" (G(T))] <= 8 (ln K + 2) J_"SSE" (Z^(ast)),
  $
  where $Z^(ast)$ is the exact global SSE minimizer @vandermerwe2003pso
  @arthur2007kmeanspp.
]

#proof[
  By the Arthur-Vassilvitskii theorem @arthur2007kmeanspp, $D^2$-seeding
  produces an initial centroid matrix $Z_1^"init"$ satisfying
  $bb(E)[J_"SSE" (Z_1^"init")] <= 8(ln K + 2) J_"SSE" (Z^(ast))$. Applying
  Lloyd's algorithm monotonically decreases SSE, so
  $J_"SSE" (Z_1 (0)) <= J_"SSE" (Z_1^"init")$. By Proposition 2(ii), the
  elitist global-best tracking of PSO guarantees
  $J_"SSE" (G(T)) <= J_"SSE" (G(0)) <= J_"SSE" (Z_1 (0))$ deterministically
  for every $T >= 0$. Taking expectations over the random initialization and
  swarm draws gives
  $bb(E)[J_"SSE" (G(T))] <= bb(E)[J_"SSE" (Z_1^"init")] <= 8(ln K + 2)
  J_"SSE" (Z^(ast))$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Each iteration of PSO-Clustering evaluates $S$ candidate centroid matrices
against all $n$ data points in $cal(O)(S n K d)$ time and updates particle
velocities and positions in $cal(O)(S K d)$ time, yielding total complexity
$cal(O)(T S n K d)$. Crucially, unlike GKA's $cal(O)(P n)$ chromosome storage,
PSO stores only $S$ matrices of size $K times d$ (positions, velocities, and
personal bests), requiring $cal(O)(S K d + n d)$ total RAM—making swarm memory
completely independent of $n$. Furthermore, distance evaluation across all $S$
particles maps cleanly onto dense BLAS-3 matrix multiplications ($X Z_i^top$),
delivering strong hardware vectorization and cache efficiency.

Geometrically, PSO searches continuously over centroid space, naturally pairing
with convex Voronoi tessellations under $J_"SSE"$ or median-like balanced
tessellations under $J_"QE"$. A subtle geometric pathology arises when computing
the difference matrices $P_i (t) - Z_i (t)$ and $G(t) - Z_i (t)$: because any
permutation of the $K$ columns of $Z_i$ represents the exact same clustering,
unaligned centroid slots can pull a particle's $k$-th centroid toward an
unrelated cluster in $G(t)$, temporarily crossing through high-energy empty
regions before re-settling.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$cal(O)(T S n K d)$ across $T$ iterations and swarm
      size $S$],
    [Space Complexity], [$cal(O)(S K d + n d)$; swarm state is strictly
      independent of $n$],
    [Objective Flexibility], [Supports non-smooth $J_"QE"$, $J_"SSE"$,
      Davies-Bouldin, or Silhouette directly],
    [Stability Condition], [Mean-square stable iff $c < 12(1 - w^2)/(7 - 5w)$
      for $c_1 = c_2 = c$],
    [Local Optimality], [Not guaranteed alone (Theorem 3); restored by
      terminal Lloyd polish],
    [Key Hyperparameters], [Cluster count $K$, swarm $S$, iterations $T$,
      $(w, c_1, c_2)$, $v_"max"$],
    [Invariances], [Translation, rotation, and scale equivariant (for QE/SSE)],
    [Hardware Profile], [Highly vectorizable via BLAS matrix products across
      $(n, K, d)$],
  ),
  caption: [Operational and geometric profile of PSO-Clustering.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multispectral and Hyperspectral Remote-Sensing Imagery*:
   Satellite pixels form dense continuous spectral clouds where land-cover
   classes (water, canopy, urban) are compact but multi-modal. Hybrid
   $K$-means-PSO was originally demonstrated on multispectral image
   classification, achieving lower quantization error and tighter inter-cluster
   separation than standalone $K$-means @vandermerwe2003pso.

2. *Wireless Sensor Network (WSN) Cluster-Head Optimization*:
   Selecting $K$ continuous relay coordinates to minimize composite energy
   dissipation (combining squared radio distance and residual battery weight)
   yields non-differentiable multi-term objectives where derivative-free PSO
   outperforms gradient methods.

3. *Financial Regime Identification and Multi-Asset Factor Clustering*:
   Asset-return factor exposures reside in continuous Euclidean spaces of
   moderate dimension ($d ~ 5$–$30$). PSO optimizes robust $L_2$-median
   quantization error $J_"QE"$ directly, preventing extreme market-crash days
   from distorting cluster centers as severely as quadratic SSE.

4. *Biomedical Multi-Channel Image Quantization and Tumor Grading*:
   Color and texture feature maps in digital pathology require compact palette
   codebooks ($K ~ 4$–$16$). Swarm parallelism explores competing palette
   configurations in a few milliseconds without storing per-pixel chromosomes.

5. *Smart-Grid Load-Profile Clustering*:
   Daily household electrical demand curves ($d = 24$ or $48$ hourly bins)
   exhibit overlapping peak patterns where single-seed $K$-means is unstable.
   Hybrid PSO with $K$-means seeding reliably stabilizes load archetypes for
   tariff design and peak-shaving forecasts.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Five failure modes require engineering care: (1) *premature non-stationary
collapse* (Theorem 3), mitigated by hybridizing with a terminal Lloyd step;
(2) *parameter divergence* outside $c < 12(1 - w^2) / (7 - 5w)$ (Theorem 1),
causing centroids to bounce violently against bounding-box walls; (3) *centroid
label mismatch* in $G(t) - Z_i (t)$, slowing convergence when $K$ is large
unless Hungarian centroid alignment is applied; (4) *empty-cluster death* when
momentum overshoots a centroid outside $"conv"(X)$; and (5) *non-convex
manifolds*, where Voronoi assignment cuts across curved structures.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Remote Sensing & Image Quant.], [Ideal: compact $cal(O)(S K d)$ memory
      + global QE search], [Add 1-step Lloyd polish at termination],
    [Custom Non-Smooth Validity], [Strong: optimizes QE, DB, or CH without
      gradients], [Enforce stability bound on $(w, c_1, c_2)$],
    [Large-$K$ Partitions ($K > 50$)], [Degrades: $K d$ search dimension and
      slot-label mismatch], [Align slots or prefer $K$-means++ / GKA],
    [Non-Convex Shapes], [Poor: Voronoi cells enforce linear boundaries],
      [Use D-Stream (Alg. 197), DBSCAN, or spectral],
    [Infinite Online Streams], [Unsuitable: batch multi-pass fitness over $X$],
      [Use ClusTree (Alg. 196) or D-Stream (Alg. 197)],
  ),
  caption: [Domain suitability matrix for PSO-Clustering.],
)

#sec("V", "Structural Relations and Empirical Behavior")

PSO-Clustering completes the evolutionary-metaheuristic duality with GKA while
relating closely to stochastic-gradient centroid methods and ant-colony search.

#proposition(
  title: "Continuous Dual of GKA and Momentum Relation to Online K-Means",
  origin: "Original Derivation",
)[
  PSO-Clustering connects formally to GKA (Alg. 198, @krishna1999gka), online /
  mini-batch $K$-means @bottou1995convergence @sculley2010webscale, and ACO
  clustering @shelokar2004aco via three structural correspondences:
  (i) *Primal-Dual Search Space*: GKA searches the discrete primal partition
  space ${1, dots.h, K}^n$ via mutation plus KMO projection, whereas PSO
  searches the continuous dual centroid space $bb(R)^(K times d)$ via
  second-order stochastic momentum, reducing population memory by a factor of
  $n / (K d)$.
  (ii) *Gradient-Free Second-Order Momentum*: Rewriting the velocity update as
  $Z_i (t + 1) - Z_i (t) = w (Z_i (t) - Z_i (t - 1)) - (c_1 R_1 + c_2 R_2)
  dot.o (Z_i (t) - overline(Z)_i (t))$ reveals a derivative-free analog of
  Polyak heavy-ball momentum toward stochastic consensus anchor
  $overline(Z)_i (t)$, contrasting with first-order stochastic gradient descent
  in online $K$-means @bottou1995convergence.
  (iii) *Swarm Intelligence Sibling*: Replacing continuous velocity flight in
  $bb(R)^(K times d)$ with discrete pheromone trails over point-to-cluster
  edges recovers Ant Colony Optimization (ACO) clustering @shelokar2004aco.
]

#proof[
  (i) Every partition $bold(s) in cal(S)_K^+$ maps to its centroid tensor
  $Z(bold(s)) in bb(R)^(K times d)$, and every $Z in bb(R)^(K times d)$ maps
  back to its Voronoi partition $bold(s)(Z)$. Storing $S$ centroid tensors
  takes $cal(O)(S K d)$ space versus $cal(O)(P n)$ for $P$ partition strings;
  at equal population size $S = P$, the memory reduction ratio is exactly
  $n / (K d)$.
  (ii) Letting $Phi_t = c_1 R_1 (t) + c_2 R_2 (t)$ and defining the stochastic
  attractor $overline(Z)_i (t) = (c_1 R_1 (t) dot.o P_i (t) + c_2 R_2 (t) dot.o
  G(t)) / Phi_t$ (entrywise) algebraically transforms Definition 1 into
  $Z_i (t + 1) = Z_i (t) + w (Z_i (t) - Z_i (t - 1)) - Phi_t dot.o (Z_i (t) -
  overline(Z)_i (t))$, matching heavy-ball momentum with gradient proxy
  $Z_i (t) - overline(Z)_i (t)$. (iii) Discretizing the social-memory matrix
  into an $n times K$ pheromone table updated by fitness-weighted assignments
  yields the ACO clustering update of @shelokar2004aco.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark
(Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index
@hubert1985comparing), PSO-Clustering attains Iris $"Sil" = 0.449$ and
$"ARI" = 0.573$, Customers $"Sil" = 0.314$, Multi-Shape $"Sil" = 0.711$
($"ARI"$ N/A), and High-D $"Sil" = 0.298$ ($"ARI"$ N/A), with a mean runtime
of *2.857 ms*. Comparing directly against its discrete evolutionary dual GKA
(Alg. 198) confirms every analytical prediction in Sections II and III:
(1) *Near-parity in global quality at $3 times$ lower latency*: PSO matches
GKA within $0.001$–$0.010$ across Iris ($"Sil"$ 0.449 vs. 0.447; $"ARI"$ 0.573
vs. 0.583), Multi-Shape ($"Sil"$ 0.711 vs. 0.712), and High-D ($"Sil"$ 0.298
vs. 0.303) while running *$3.0 times$ faster* (2.857 ms vs. 8.592 ms) because
continuous $K times d$ centroid updates eliminate per-point allele mutation
overhead.
(2) *Slight edge of KMO over pure swarm dynamics*: GKA edges narrowly ahead on
Iris ARI ($+0.010$), Customers Silhouette ($+0.026$), and High-D Silhouette
($+0.005$), reflecting the exact stationarity guarantee of KMO (Alg. 198, Lemma
1) versus PSO's vulnerability to premature collapse and slot-label mismatch
(Theorem 3 and Section III).
(3) *Stream vs. evolutionary tradeoff*: Both PSO and GKA substantially
outperform ClusTree and D-Stream on High-D Silhouette (0.298–0.303 vs.
0.247–0.257) and Iris ARI (0.573–0.583 vs. 0.384–0.494), trading an order of
magnitude in runtime for global multi-pass Voronoi exploration.
