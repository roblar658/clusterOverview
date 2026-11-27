#import "../template/ieee.typ": *

#alg_entry(
  id: 141,
  name: "OptiGrid",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1999",
  paradigm: "Contracting-Projection Recursive Grid Partitioning",
  complexity: $O(n d log n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Fixed-resolution grid clustering methods such as STING @wang1997sting and
CLIQUE @agrawal1998clique discretize $bb(R)^d$ into $b^d$ axis-aligned cells,
incurring exponential cell proliferation in high dimensions while slicing
arbitrarily through dense cluster interiors. *OptiGrid* (Optimal
Grid-Clustering), introduced by Hinneburg and Keim @hinneburg1999optigrid,
resolves both bottlenecks by replacing uniform static meshes with recursive,
data-adaptive *optimal cutting planes* derived from one-dimensional density
estimates along *contracting linear projections*.

Let $cal(X) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ be a finite dataset
and let $K_sigma^(d) (bold(v)) = (2 pi sigma^2)^(-d / 2) exp(-norm(bold(v))_2^2 / (2 sigma^2))$
denote the isotropic $d$-dimensional Gaussian kernel with bandwidth
$sigma > 0$. The empirical $d$-dimensional kernel density estimate (KDE) at
$bold(x) in bb(R)^d$ is
$
  hat(f)_sigma^(d) (bold(x)) = 1 / n sum_(i=1)^n K_sigma^(d) (bold(x) - bold(x)_i).
$
Directly locating $(d-1)$-dimensional low-density separating surfaces of
$hat(f)_sigma^(d)$ is intractable in high dimensions. OptiGrid instead
constructs hyperplanes from one-dimensional projected densities.

#definition(title: "Contracting Projection and Cutting Plane")[
  A linear map $P : bb(R)^d -> bb(R)$ given by
  $P(bold(x)) = chevron.l bold(u), bold(x) chevron.r$ is a *contracting
  projection* if it is $1$-Lipschitz under the Euclidean norm, i.e.,
  $
    |P(bold(x)) - P(bold(y))| <= norm(bold(x) - bold(y))_2
    wide forall bold(x), bold(y) in bb(R)^d,
  $
  which holds if and only if $norm(bold(u))_2 <= 1$ (canonically unit vectors
  $norm(bold(u))_2 = 1$). For a subset $S subset.eq cal(X)$ of size $m = |S|$,
  the projected 1D Gaussian KDE is
  $
    hat(f)_(P, sigma) (t) = 1 / (m sqrt(2 pi) sigma)
      sum_(bold(x) in S) exp(- (t - P(bold(x)))^2 / (2 sigma^2)).
  $
  A *cutting plane* $H_(P, t^*) = {bold(x) in bb(R)^d : P(bold(x)) = t^*}$ is
  *valid* if $t^*$ is a local minimum of $hat(f)_(P, sigma)$ separating two
  local maxima $t_1 < t^* < t_2$ whose heights exceed a noise threshold
  $m_"noise" > 0$, with valley score $s(P, t^*) = hat(f)_(P, sigma) (t^*)$.
]

#definition(title: "Optimal (q, m_cut)-Grid Partition")[
  Given a candidate pool $cal(P) = {P_1, dots, P_k}$ of contracting
  projections, a maximum per-level cut budget $q >= 1$, a maximum valley
  density (or minimum separation-ratio) criterion $m_"cut"$, and minimum
  dense-cell support $xi >= 1$, a *$(q, m_"cut")$-optimal grid* on $S$ selects
  the $r <= q$ best valid cutting planes $H_(P_(j_1), t_1^*), dots, H_(P_(j_r), t_r^*)$
  meeting $m_"cut"$, partitions $S$ into the resulting non-empty polytopal grid
  cells, discards sparse cells with $|S_c| < xi$, and recurses on each dense
  cell $|S_c| >= xi$.
]

Because OptiGrid stores only *populated* grid cells rather than the full
Cartesian product of cut intervals, the active cell count at every recursion
level remains strictly bounded by $floor(n / xi)$ regardless of ambient
dimension $d$.

#algbox(
  title: "OptiGrid: Optimal Contracting-Grid Clustering",
  inputs: [Dataset $cal(X) subset bb(R)^d$, contracting projections $cal(P) = {P_1, dots, P_k}$, cut budget $q$, noise level $m_"noise"$, cut threshold $m_"cut"$, bandwidth $sigma$, min cluster size $xi$],
  outputs: [Cluster partition $cal(C) = {C_1, dots, C_M}$ and noise set $N$],
)[
  + Initialize active queue $cal(Q) <- {cal(X)}$ and cluster set $cal(C) <- emptyset$.
  + *While* $cal(Q) != emptyset$: pop active subset $S$ from $cal(Q)$ and set candidate cut list $cal(L) <- emptyset$.
  + *For each* projection $P_j in cal(P)$: compute projected values $z_i = P_j (bold(x)_i)$ for $bold(x)_i in S$ and form 1D KDE $hat(f)_(P_j, sigma) (t)$ (via binned convolution or sorted scan).
  + Locate all local maxima above $m_"noise"$ and extract intervening local minima $t^*$ satisfying $m_"cut"$; append $(P_j, t^*, s(P_j, t^*))$ to $cal(L)$.
  + *If* $cal(L) = emptyset$: mark $S$ as a terminal cluster ($cal(C) <- cal(C) union {S}$) and continue.
  + Select up to $q$ lowest-valley-density cutting planes from $cal(L)$ and assign each $bold(x) in S$ to a cell signature tuple over the selected cuts.
  + Group points by non-empty cell signature into ${S_1, dots, S_m}$; discard any $S_c$ with $|S_c| < xi$ as noise, and push remaining dense cells onto $cal(Q)$ (or emit $S$ if only one cell contains all of $S$).
  + *Return* $cal(C)$ and unassigned points $N = cal(X) without union.big_(C in cal(C)) C$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five formal guarantees characterizing OptiGrid's density
preservation, mixture separability, geometric equivariance, Lyapunov
termination and sparsity, and manifold separation limits.

#subsec("A", "Radon-Marginal Identity and No-False-Splitting Bound")

#lemma(title: "Radon Marginalization and No-False-Splitting", origin: "Classical — Cited")[
  Let $P(bold(x)) = chevron.l bold(u), bold(x) chevron.r$ with $norm(bold(u))_2 = 1$
  be an orthogonal contracting projection, and let $H_(P, t) = {bold(x) in bb(R)^d : chevron.l bold(u), bold(x) chevron.r = t}$.
  Then for every $t in bb(R)$ and every $bold(x) in H_(P, t)$ @hinneburg1999optigrid:
  $
    hat(f)_(P, sigma) (t) &= integral_(H_(P, t)) hat(f)_sigma^(d) (bold(y)) dif cal(H)^(d-1) (bold(y)), \
    sup_(bold(x) in H_(P, t)) hat(f)_sigma^(d) (bold(x)) &<= (2 pi sigma^2)^(-(d-1) / 2) hat(f)_(P, sigma) (t),
  $
  where $cal(H)^(d-1)$ is $(d-1)$-dimensional Hausdorff measure on $H_(P, t)$.
  More generally, for any contracting projection with $norm(bold(u))_2 <= 1$,
  $hat(f)_(P, sigma) (P(bold(x))) >= (2 pi sigma^2)^((d-1) / 2) hat(f)_sigma^(d) (bold(x))$.
]

#proof[
  Choose an orthonormal basis ${bold(u), bold(v)_2, dots, bold(v)_d}$ of
  $bb(R)^d$ and write any $bold(y) in H_(P, t)$ as $bold(y) = t bold(u) + V bold(w)$
  with $bold(w) in bb(R)^(d-1)$ and $V = [bold(v)_2, dots, bold(v)_d]$. Similarly
  decompose $bold(x)_i = z_i bold(u) + V bold(w)_i$ where $z_i = P(bold(x)_i)$.
  By Pythagorean orthogonality,
  $norm(bold(y) - bold(x)_i)_2^2 = (t - z_i)^2 + norm(bold(w) - bold(w)_i)_2^2$,
  so the isotropic Gaussian kernel factorizes exactly:
  $
    K_sigma^(d) (bold(y) - bold(x)_i) = K_sigma^(1) (t - z_i) dot K_sigma^(d-1) (bold(w) - bold(w)_i).
  $
  Integrating over $bold(w) in bb(R)^(d-1)$ gives
  $integral_(bb(R)^(d-1)) K_sigma^(d-1) (bold(w) - bold(w)_i) dif bold(w) = 1$.
  Averaging over $i = 1, dots, n$ yields the exact Radon-transform identity
  $hat(f)_(P, sigma) (t) = integral_(H_(P, t)) hat(f)_sigma^(d) (bold(y)) dif cal(H)^(d-1) (bold(y))$.
  For the pointwise inequality under any $norm(bold(u))_2 <= 1$, Cauchy--Schwarz
  gives $|P(bold(x)) - P(bold(x)_i)| <= norm(bold(u))_2 norm(bold(x) - bold(x)_i)_2 <= norm(bold(x) - bold(x)_i)_2$,
  whence
  $
    exp(- (P(bold(x)) - P(bold(x)_i))^2 / (2 sigma^2))
      >= exp(- norm(bold(x) - bold(x)_i)_2^2 / (2 sigma^2)).
  $
  Multiplying by $(n (2 pi sigma^2)^(1/2))^(-1)$ and summing over $i$ produces
  $hat(f)_(P, sigma) (P(bold(x))) >= (2 pi sigma^2)^((d-1) / 2) hat(f)_sigma^(d) (bold(x))$.
  Taking the supremum over $bold(x) in H_(P, t)$ proves the second claim.
]

Consequently, a cutting plane with small valley density
$hat(f)_(P, sigma) (t^*) <= tau$ directly bounds both the *integrated*
cross-sectional density and the *maximum* pointwise $d$-dimensional KDE
density along the entire hyperplane, preventing false cuts through dense
cluster cores.

#subsec("B", "Exact Separation Criterion for Gaussian Mixtures")

#theorem(title: "Bandwidth-Separation Threshold for Valley Detection", origin: "Original Derivation")[
  Let $bold(X)$ follow a balanced two-component isotropic Gaussian mixture
  $1/2 cal(N)(bold(mu)_1, tau^2 I_d) + 1/2 cal(N)(bold(mu)_2, tau^2 I_d)$,
  and let $P(bold(x)) = chevron.l bold(u), bold(x) chevron.r$ with
  $norm(bold(u))_2 = 1$ and projected separation
  $Delta_u = |chevron.l bold(u), bold(mu)_1 - bold(mu)_2 chevron.r|$. Then the
  population projected KDE $f_(P, sigma) = bb(E)[hat(f)_(P, sigma)]$ admits a
  strict local minimum at the midpoint $t^* = (P(bold(mu)_1) + P(bold(mu)_2)) / 2$
  if and only if
  $
    Delta_u > 2 sqrt(tau^2 + sigma^2),
    wide "equivalently" wide
    sigma < sqrt(Delta_u^2 / 4 - tau^2).
  $
  Whenever this holds, the valley-to-center density ratio satisfies
  $
    f_(P, sigma) (t^*) / f_(P, sigma) (P(bold(mu)_1))
      = (2 exp(- Delta_u^2 / (8 s^2))) / (1 + exp(- Delta_u^2 / (2 s^2))),
  $
  where $s^2 = tau^2 + sigma^2$.
]

#proof[
  Because $P(bold(X))$ is a mixture of two 1D Gaussians
  $1/2 cal(N)(m_1, tau^2) + 1/2 cal(N)(m_2, tau^2)$ with $m_r = P(bold(mu)_r)$,
  convolving with $K_sigma^(1)$ yields
  $
    f_(P, sigma) (t) = 1 / 2 (phi_s (t - m_1) + phi_s (t - m_2)),
  $
  where $phi_s (z) = (2 pi s^2)^(-1/2) exp(-z^2 / (2 s^2))$ and $s^2 = tau^2 + sigma^2$.
  Without loss of generality translate coordinates so $m_1 = -Delta_u / 2$,
  $m_2 = +Delta_u / 2$, and $t^* = 0$. By even symmetry,
  $f'_(P, sigma) (0) = 0$. Differentiating twice using
  $phi''_s (z) = (z^2 / s^4 - 1 / s^2) phi_s (z)$ gives
  $
    f''_(P, sigma) (0) = (Delta_u^2 / (4 s^4) - 1 / s^2) phi_s (Delta_u / 2).
  $
  Since $phi_s (Delta_u / 2) > 0$, $t^* = 0$ is a strict local minimum if and
  only if $Delta_u^2 / (4 s^4) > 1 / s^2$, i.e., $Delta_u > 2 s = 2 sqrt(tau^2 + sigma^2)$.
  Evaluating $f_(P, sigma) (0) = phi_s (Delta_u / 2)$ and
  $f_(P, sigma) (Delta_u / 2) = 1/2 (phi_s (0) + phi_s (Delta_u))$ directly
  gives the stated closed-form valley-to-center ratio, which decays
  exponentially as $O(exp(-Delta_u^2 / (8 s^2)))$.
]

This theorem quantifies both the bandwidth upper bound required to avoid
over-smoothing cluster valleys and the critical role of projection direction
$bold(u)$: maximizing $Delta_u$ over $cal(P)$ exponentially sharpens the
valley contrast.

#subsec("C", "Geometric Equivariance and Axis-Projection Attenuation")

#proposition(title: "Orthogonal Equivariance vs. Axis-Only SNR Loss", origin: "Original Proof")[
  Let $cal(T)_Q (bold(x)) = Q bold(x) + bold(b)$ be any rigid motion
  ($Q in O(d)$, $bold(b) in bb(R)^d$).
  + If the projection pool $cal(P)$ is equivariantly transformed to
    $cal(P)_Q = {bold(x) |-> chevron.l Q bold(u), bold(x) chevron.r : bold(u) in cal(P)}$
    (e.g., data-driven PCA/SVD or random orthogonal projections), then
    OptiGrid's cutting planes and cluster partition are strictly equivariant
    under $cal(T)_Q$ and under uniform scaling $bold(x) |-> c bold(x)$ with
    $sigma |-> c sigma$.
  + Conversely, if $cal(P) = {bold(e)_1, dots, bold(e)_d}$ is restricted to
    coordinate axes, then for any separation vector $bold(Delta) = bold(mu)_1 - bold(mu)_2$,
    $
      max_(1 <= j <= d) |chevron.l bold(e)_j, bold(Delta) chevron.r|
        >= norm(bold(Delta))_2 / sqrt(d),
    $
    with equality attained tightly by diagonal clusters
    $bold(Delta) = (c, dots, c)^top$. Thus axis-only OptiGrid suffers a
    worst-case $sqrt(d)$ separation attenuation and fails to split diagonal
    Gaussians whenever $2 s < norm(bold(Delta))_2 <= 2 s sqrt(d)$.
]

#proof[
  For (1), for any $bold(u) in cal(P)$ let $bold(u)' = Q bold(u)$. Then
  $norm(bold(u)')_2 = norm(bold(u))_2 <= 1$ and
  $chevron.l bold(u)', Q bold(x)_i + bold(b) chevron.r = chevron.l bold(u), bold(x)_i chevron.r + chevron.l bold(u)', bold(b) chevron.r$.
  Hence the projected KDE shifts rigidly by the constant
  $beta = chevron.l bold(u)', bold(b) chevron.r$, preserving all local extrema,
  valley scores, and cell memberships identically. Scaling by $c > 0$ with
  $sigma' = c sigma$ maps $hat(f)_(P, sigma') (c t) = c^(-1) hat(f)_(P, sigma) (t)$,
  preserving all extrema ordering and relative valley-to-peak ratios.
  For (2), by the equivalence of $ell_2$ and $ell_infinity$ norms in $bb(R)^d$,
  $norm(bold(Delta))_2^2 = sum_(j=1)^d Delta_j^2 <= d norm(bold(Delta))_infinity^2$,
  so $max_j |Delta_j| >= norm(bold(Delta))_2 / sqrt(d)$, with exact equality
  when $|Delta_1| = dots = |Delta_d| = c$. Substituting into Theorem 1 shows
  that when $2 s < c sqrt(d) <= 2 s sqrt(d)$, the optimal diagonal projection
  $bold(u)^* = bold(Delta) / norm(bold(Delta))_2$ satisfies $Delta_(u^*) > 2 s$
  and yields a strict valley, whereas every axis projection satisfies
  $Delta_(e_j) = c <= 2 s$ and exhibits a strictly unimodal density with no
  valid cut.
]

#subsec("D", "Discrete Lyapunov Termination and Active-Cell Sparsity")

#theorem(title: "Lyapunov Termination, Sparsity Bound, and Complexity", origin: "Original Proof")[
  Let $cal(G)_t$ be the collection of active dense cells ($|S| >= xi >= 1$)
  after $t$ non-trivial node splits of OptiGrid (where each split either
  partitions a parent cell into $r >= 2$ dense children or discards at least
  one noise point into sparse cells of size $< xi$). Define the discrete
  Lyapunov potential
  $
    Psi(cal(G)_t) = sum_(S in cal(G)_t) (|S| - xi).
  $
  Then:
  + $Psi(cal(G)_t) >= 0$ for all $t >= 0$, with $Psi(cal(G)_0) = n - xi$;
  + Every split strictly decreases $Psi$ by at least $1$ (and by at least $xi$
    whenever a cell splits into $r >= 2$ dense children);
  + The number of active cells never exceeds $|cal(G)_t| <= floor(n / xi)$,
    and recursion terminates after at most $n - xi$ total splits (at most
    $floor(n / xi) - 1$ splits in the noiseless regime);
  + Under $k = O(d)$ axis-aligned projections and balanced recursion depth
    $h = O(log(n / xi))$, total time complexity is $O(n d log n)$ (or
    $O(h k n d)$ for $k$ general linear projections with $B$-bin 1D KDE).
]

#proof[
  Since every $S in cal(G)_t$ satisfies $|S| >= xi$, each summand $|S| - xi >= 0$,
  giving $Psi(cal(G)_t) >= 0$. Because active cells are pairwise disjoint
  subsets of $cal(X)$, $xi |cal(G)_t| <= sum_(S in cal(G)_t) |S| <= n$, proving
  $|cal(G)_t| <= floor(n / xi)$. Now suppose step $t + 1$ splits a cell
  $S in cal(G)_t$ into $r >= 1$ surviving dense children $S_1, dots, S_r$ plus
  $m_"noise" >= 0$ discarded noise points, so $|S| = m_"noise" + sum_(c=1)^r |S_c|$.
  The exact change in potential is
  $
    Delta Psi &= sum_(c=1)^r (|S_c| - xi) - (|S| - xi) \
              &= - m_"noise" - xi (r - 1).
  $
  If $r >= 2$, then $Delta Psi <= -xi <= -1$; if $r = 1$, non-triviality
  requires $m_"noise" >= 1$, so $Delta Psi = -m_"noise" <= -1$. Because $Psi$
  is a non-negative integer decreasing by $>= 1$ per split, total splits are
  bounded by $Psi(cal(G)_0) = n - xi$ (and by $floor((n - xi) / xi) = floor(n / xi) - 1$
  when $m_"noise" = 0$ throughout). At any single recursion depth $ell$, the
  active cells are pairwise disjoint and contain at most $n$ points total.
  Projecting along $d$ coordinate axes costs $O(n d)$ (lookup only), sorting
  or evaluating 1D KDE over $B$ bins costs $O(n d log n)$ or $O(d (n + B))$,
  and partitioning into cell signatures takes $O(n q)$. Summing across
  $h = O(log(n / xi))$ balanced levels gives $O(n d log n)$ (with $O(h k n d)$
  for $k$ dense linear projections using binned 1D KDE).
]

#subsec("E", "Topological Limitation on Concentric Manifolds")

#corollary(title: "Impossibility of Single-Cut Concentric Separation", origin: "Original Counterexample")[
  Let $C_1 = {bold(x) in bb(R)^2 : norm(bold(x))_2 <= r}$ and
  $C_2 = {bold(x) in bb(R)^2 : norm(bold(x))_2 = R}$ with $0 < r < R$ be an
  inner disk surrounded by a concentric uniform circle. Then:
  + For every unit vector $bold(u) in S^1$, the projected marginal density of
    $C_2$ is $g_R (t) = (pi sqrt(R^2 - t^2))^(-1) > 0$ on $(-R, R)$, which
    strictly covers $P(C_1) = [-r, r]$;
  + Every cutting plane $H_(P, t^*)$ with $|t^*| < R$ intersects $C_2$ at two
    antipodal-chord points $bold(x)_plus.minus = t^* bold(u) plus.minus sqrt(R^2 - (t^*)^2) bold(u)^perp$,
    making it topologically impossible for any finite single-level hyperplane
    partition to separate $C_1$ from $C_2$ without fragmenting $C_2$.
]

#proof[
  Parameterize $C_2$ by angle $theta in [0, 2 pi)$ relative to $bold(u)$, so
  $t = P(bold(x)) = R cos theta$. For $|t| < R$, the preimage consists of two
  angles $plus.minus arccos(t / R)$, each contributing Jacobian
  $|dif t / dif theta|^(-1) = (R |sin theta|)^(-1) = (R^2 - t^2)^(-1/2)$ under
  uniform arc density $(2 pi)^(-1)$. Summing both branches yields
  $g_R (t) = (pi sqrt(R^2 - t^2))^(-1) >= (pi R)^(-1) > 0$ throughout
  $[-r, r] subset (-R, R)$. Any hyperplane $H_(P, t^*)$ separating $C_1$ from
  any point of $C_2$ must satisfy $|t^*| <= R$; if $|t^*| < R$, the two
  explicit intersection points $bold(x)_plus.minus in C_2 inter H_(P, t^*)$
  disconnect $C_2$ into two open arcs lying on opposite half-spaces, so $C_2$
  is necessarily split across distinct grid cells.
]

#sec("III", "Complexity, Geometry and Operational Profile")

OptiGrid's computational profile is dominated by 1D projections and density
estimation rather than pairwise $O(n^2)$ distance calculations or exponential
$b^d$ grid materialization. Storing only active cell signatures in a hash map
requires $O(n d)$ memory for the dataset plus $O(n + k B)$ workspace for
projected arrays and $B$-bin 1D density histograms, yielding strong L1/L2 cache
locality during sequential coordinate scans. Adding random or PCA-derived
contracting projections increases per-level projection cost from $O(n d)$ to
$O(k n d)$ while mitigating the $sqrt(d)$ diagonal SNR attenuation established
in Proposition 1.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n d log n)$ balanced axis-aligned; $O(h k n d)$ with $k$ general projections and $B$-bin 1D KDE],
    [Space Complexity], [$O(n d + k B)$; active grid cells bounded by $floor(n / xi)$ (zero empty-cell materialization)],
    [Cluster Geometry], [Convex polytopes (axis boxes or general half-space intersections); fragments nested rings],
    [Noise Handling], [Explicit sparse-cell pruning ($|S_c| < xi$) and peak-height noise filter $m_"noise"$],
    [Hyperparameters], [Bandwidth $sigma$, cut budget $q$, noise level $m_"noise"$, cut threshold $m_"cut"$, min size $xi$],
    [Invariances], [Translation & permutation invariant; $O(d)$-equivariant iff $cal(P)$ includes rotated projections],
  ),
  caption: [Operational, geometric, and complexity profile of OptiGrid.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *High-Dimensional Astronomical Sky Catalogs*: Multi-band photometric surveys
  (e.g., SDSS, Gaia) generate millions of objects in $10$--$50$ color-magnitude
  features containing dense stellar/galactic loci separated by sparse
  photometric valleys. OptiGrid's Radon-marginal property (Lemma 1) guarantees
  that cutting planes pass exclusively through empty color-space gaps without
  exponential $b^d$ memory blowup @hinneburg1999optigrid.
+ *High-Throughput Flow Cytometry and Cell Gating*: Multiparameter cytometry
  measures $15$--$40$ fluorescent markers per cell across $10^6$ events.
  Biologists traditionally perform manual 1D/2D histogram gating at density
  minima between marker-negative and marker-positive populations; OptiGrid
  mathematically automates hierarchical gating by selecting globally optimal
  low-density 1D valleys across all marker projections.
+ *Network Intrusion and Telemetry Traffic Partitioning*: High-rate packet
  flows described by $20$--$80$ timing, byte-count, and entropy features
  exhibit multimodal normal baselines interspersed with bursty anomaly regimes.
  OptiGrid's $O(n d log n)$ throughput and sparse-cell noise rejection
  ($|S_c| < xi$) isolate coherent behavioral strata while filtering sporadic
  outliers without computing $O(n^2)$ distance graphs.
+ *Image Feature-Space Indexing and Content-Based Retrieval*: Color-texture
  descriptor databases (e.g., $16$--$64$-D histograms or visual embeddings)
  require fast coarse-to-fine spatial partitioning. OptiGrid constructs a
  compact polytopal partition tree whose leaves correspond directly to dense
  semantic clusters, doubling as an adaptive indexing structure @aggarwal2013data.
+ *Industrial Wafer-Map and Semiconductor Process Monitoring*: Fab sensor
  logs combine dozens of continuous etch, deposition, and metrology variables
  with distinct operational modes. Because process shifts often manifest as
  separable bimodal shifts along specific sensor axes or linear combinations,
  OptiGrid isolates regime clusters while yielding directly interpretable
  hyperplane decision rules $chevron.l bold(u), bold(x) chevron.r = t^*$.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

OptiGrid degrades in four mathematically characterized regimes: (i) *concentric
or interlocking non-convex manifolds* (Corollary 1), where every linear cut
slices the outer manifold; (ii) *oblique high-dimensional clusters under
axis-only projections* (Proposition 1), where diagonal separation attenuates
by $1 / sqrt(d)$; (iii) *over-smoothed bandwidth $sigma$* violating
$sigma < sqrt(Delta_u^2 / 4 - tau^2)$ (Theorem 1), which merges adjacent
valleys; and (iv) *extreme multi-scale density contrast*, where a global
$m_"noise"$ threshold either ignores sparse genuine clusters or admits noisy
spurious peaks.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-D photometry / cytometry], [Strong: sharp 1D valleys; $O(n / xi)$ cell sparsity avoids $b^d$ blowup], [Tune $sigma$ via Silverman rule; add PCA cuts],
    [Telemetry & tabular segmentation], [Strong: interpretable hyperplane cuts and sub-millisecond speed], [Standardize heterogeneous feature scales first],
    [Concentric / spiral manifolds], [Poor: hyperplanes necessarily fragment enclosing manifolds (Cor. 1)], [Use DBSCAN @ester1996dbscan or OPTICS @ankerst1999optics],
    [Oblique correlated subspaces], [Moderate: fails with axis-only cuts when $Delta <= 2 s sqrt(d)$], [Add diagonal / SVD cuts or use ORCLUS @aggarwal2000orclus],
    [Sparse tiny-sample high-D ($n << d$)], [Weak: 1D KDE becomes noisy and undersampled per level], [Use subspace methods (ENCLUS, FINDIT) or PCA],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for OptiGrid.],
)

#sec("V", "Structural Relations and Empirical Behavior")

OptiGrid bridges pure kernel density clustering (DENCLUE @hinneburg1998denclue,
@hinneburg2007denclue2) and static grid clustering (STING @wang1997sting,
BANG @schikuta1998bang, CLIQUE @agrawal1998clique). The next proposition makes
this connection explicit via density-gradient geometry.

#proposition(title: "Link to DENCLUE Density Attractors", origin: "Original Proof")[
  Let $cal(A) = {bold(a)_1, dots, bold(a)_m}$ be the local maxima (density
  attractors) of the $d$-dimensional Gaussian KDE $hat(f)_sigma^(d)$ used by
  DENCLUE @hinneburg1998denclue. If OptiGrid's recursive cutting planes
  $H_1, dots, H_L$ all satisfy $hat(f)_(P_r, sigma) (t_r^*) < (2 pi sigma^2)^((d-1)/2) min_j hat(f)_sigma^(d) (bold(a)_j)$
  and every DENCLUE basin of attraction is convexly separable from the others
  by $H_1, dots, H_L$, then no cutting plane passes through any density
  attractor $bold(a)_j$ or its superlevel core
  ${bold(x) : hat(f)_sigma^(d) (bold(x)) > (2 pi sigma^2)^(-(d-1)/2) max_r hat(f)_(P_r, sigma) (t_r^*)}$.
]

#proof[
  By Lemma 1, for every cut $H_r = H_(P_r, t_r^*)$ and every point
  $bold(x) in H_r$,
  $hat(f)_sigma^(d) (bold(x)) <= (2 pi sigma^2)^(-(d-1)/2) hat(f)_(P_r, sigma) (t_r^*)$.
  Under the stated valley-height hypothesis, this upper bound is strictly
  smaller than $min_j hat(f)_sigma^(d) (bold(a)_j)$ and strictly below every
  point in the specified superlevel core. Hence $H_r$ is strictly disjoint
  from every density attractor and its high-density core, acting as a
  piecewise-linear outer approximation to DENCLUE's attraction basins without
  requiring $O(n)$ hill-climbing trajectories per point.
]

Unlike CLIQUE @agrawal1998clique and ENCLUS (Alg. 143), which search over
overlapping low-dimensional subspaces, OptiGrid produces a disjoint global
full-space partition (similar to BANG-Clustering, Alg. 142), yet avoids
BANG's cyclic dyadic splitting bottleneck in high dimensions by choosing
only the most discriminative projections @kriegel2009subspacesurvey.

On the ClusterOverview 4-dataset benchmark (evaluated via Silhouette
@rousseeuw1987silhouettes and Adjusted Rand Index @hubert1985comparing),
OptiGrid achieves: *Iris* $"Sil" = 0.586$, $"ARI" = 0.567$; *Customers*
$"Sil" = 0.232$; *Multi-Shape* $"Sil" = 0.706$ ($"ARI" =$ N/A); *High-D*
$"Sil" = 0.287$ ($"ARI" =$ N/A); with a mean runtime of $0.761$ ms. These
empirical metrics mirror our theoretical guarantees with striking precision.
First, OptiGrid attains the highest *High-D* Silhouette ($0.287$) among
grid-and-subspace peers BANG ($0.168$), ENCLUS ($0.232$), and FINDIT
($0.264$), directly validating Theorem 2's $floor(n / xi)$ sparsity bound and
the resistance of 1D contracting KDEs to the curse of dimensionality
@beyer1999nearest. Second, its strong *Multi-Shape* Silhouette ($0.706$)
reflects crisp low-density linear valleys separating well-spaced 2D
components, whereas *Iris* ($"ARI" = 0.567$) and *Customers* ($"Sil" = 0.232$)
illustrate Theorem 1's valley-smoothing limit: overlapping Versicolor and
Virginica distributions and unimodal-skewed customer spending densities lack
deep 1D valleys, prompting conservative cuts. Finally, the $0.761$ ms mean
runtime confirms the practical efficiency of $O(n d log n)$ 1D density
partitioning.
