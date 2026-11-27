#import "../template/ieee.typ": *

#alg_entry(
  id: 144,
  name: "FINDIT",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2004",
  paradigm: "Medoid-Sampled Dimension-Oriented Distance Subspace Clustering",
  complexity: $O(n k d + |S| |M| d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Projected clustering algorithms such as PROCLUS @aggarwal1999proclus partition
high-dimensional data into disjoint clusters, each associated with its own
subset of correlated coordinates, thereby overcoming the exponential subspace
explosion of grid miners like CLIQUE @agrawal1998clique and ENCLUS
@cheng1999enclus. However, PROCLUS requires two notoriously difficult user
inputs---the exact cluster count $k$ and average subspace dimensionality
$l$---and evaluates candidate medoids using $L_1$ distances that suffer from
high-dimensional distance concentration @beyer1999nearest. *FINDIT* (Fast and
INtelligent subspace clustering algorithm using DImension voTing), introduced
by Woo, Lee, Kim, and Lee @woo2004findit, resolves both limitations via two
mathematical innovations: the *Dimension-Oriented Distance (DOD)* and
*Dimension Voting* over Chernoff-sized random samples.

Let $cal(X) = {bold(x)_1, dots, bold(x)_n} subset bb(R)^d$ be a standardized
dataset and let $epsilon > 0$ be a per-coordinate proximity threshold (set
directly or adapted from coordinate dispersion).

#definition(title: "Dimension-Oriented Distance (DOD)")[
  For any two points $bold(x), bold(y) in bb(R)^d$ and tolerance $epsilon > 0$,
  define their *correlated dimension set*
  $
    D_epsilon (bold(x), bold(y)) = {j in {1, dots, d} : |x_j - y_j| <= epsilon}.
  $
  The *Dimension-Oriented Distance (DOD)* @woo2004findit between $bold(x)$
  and $bold(y)$ counts the number of disagreeing coordinates:
  $
    d_"DOD"^(epsilon) (bold(x), bold(y)) = d - |D_epsilon (bold(x), bold(y))|
      = sum_(j=1)^d bb(I)(|x_j - y_j| > epsilon).
  $
]

Unlike Euclidean or Manhattan distance---where a single noisy irrelevant
coordinate can contribute an unbounded penalty---DOD bounds the contribution
of every irrelevant dimension strictly by $1$, preserving strong proximity
between points sharing a correlated subspace $D_r subset.eq {1, dots, d}$.

#definition(title: "Dimension Voting and Subspace Assignment")[
  Given a candidate medoid $bold(m)_r$ and its DOD-neighborhood
  $N_r subset.eq S$ in a random sample $S subset cal(X)$, each neighbor
  $bold(x) in N_r$ casts a coordinate vote
  $v_j (bold(x), bold(m)_r) = bb(I)(|x_j - m_(r, j)| <= epsilon)$ for each
  $j in {1, dots, d}$. The *voted correlated subspace* of cluster $r$ at
  support threshold $tau_"vote" in (0, 1)$ is
  $
    hat(D)_r = {j in {1, dots, d} :
      1 / |N_r| sum_(bold(x) in N_r) v_j (bold(x), bold(m)_r) >= tau_"vote"}.
  $
  In the final assignment phase, each point $bold(x) in cal(X)$ is assigned to
  the cluster $r^*$ minimizing the subspace-normalized Manhattan distance
  $
    r^* (bold(x)) = arg min_(r in {1, dots, k})
      1 / |hat(D)_r| sum_(j in hat(D)_r) |x_j - m_(r, j)|,
  $
  or labeled as an outlier if its minimum distance exceeds the cluster's
  corridor bound.
]

FINDIT structures execution into three clean phases: (i) *Sampling*, drawing
medoid candidate set $M$ and evaluation sample $S$ of $n$-independent sizes;
(ii) *Cluster Forming*, iteratively expanding and merging DOD neighborhoods
around medoids in $M$, tuning $epsilon$ to maximize cluster quality, and
extracting each cluster's subspace $hat(D)_r$ via dimension voting; and
(iii) *Data Assignment*, classifying all $n$ points in a single linear scan.

#algbox(
  title: "FINDIT: Subspace Clustering via DOD and Dimension Voting",
  inputs: [Standardized dataset $cal(X) subset bb(R)^d$, min cluster fraction $alpha_"min"$, minimum subspace size $l_"min"$, tolerance grid $cal(E)$, voting threshold $tau_"vote"$],
  outputs: [Projected clusters ${(C_r, bold(m)_r, hat(D)_r)}_(r=1)^k$ and outlier set $N$],
)[
  + *Sampling Phase*: Draw independent uniform random samples $S, M subset cal(X)$ sized via Chernoff bounds for $alpha_"min"$ (typically $|S| >= |M|$, independent of $n$).
  + *For each* candidate threshold $epsilon in cal(E)$: compute pairwise DOD matrix $d_"DOD"^(epsilon) (bold(m), bold(x))$ for $bold(m) in M, bold(x) in S$.
  + Form initial neighborhoods $N_epsilon (bold(m)) = {bold(x) in S : d_"DOD"^(epsilon) (bold(m), bold(x)) <= d - l_"min"}$ and prune sparse medoids with $|N_epsilon (bold(m))| < c_"min"$.
  + Greedily select highest-support medoids, merging overlapping neighborhoods whose medoids agree on $>= l_"min"$ dimensions (small mutual DOD) into candidate clusters ${N_1, dots, N_k}$.
  + *Dimension Voting*: For each $N_r$, tally coordinate votes $hat(p)_(r, j) = |N_r|^(-1) sum_(bold(x) in N_r) bb(I)(|x_j - m_(r, j)| <= epsilon)$ and select subspace $hat(D)_r = {j : hat(p)_(r, j) >= tau_"vote"}$.
  + Retain the $epsilon^*$ configuration maximizing intra-subspace compactness and inter-cluster separation.
  + *Data Assignment Phase*: In one $O(n k d)$ scan, assign each $bold(x) in cal(X)$ to $r^* (bold(x)) = arg min_r |hat(D)_r|^(-1) sum_(j in hat(D)_r) |x_j - m_(r, j)|$ if within cluster radius, else to outlier set $N$.
  + *Return* ${(C_r, bold(m)_r, hat(D)_r)}_(r=1)^k$ and $N$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five rigorous results covering DOD's relaxed metric geometry,
Bernstein high-dimensional contrast preservation, exact Chernoff sample-size
bounds, finite-sample subspace recovery under majority voting, and coordinate
breakdown analysis.

#subsec("A", "Metric Axiom Audit and Scale-Doubled Triangle Inequality")

#proposition(title: "Semi-Metric Structure and 2-Epsilon Triangle Inequality", origin: "Original Proof")[
  For any $epsilon > 0$, the Dimension-Oriented Distance $d_"DOD"^(epsilon)$
  on $bb(R)^d$ is non-negative, reflexive, and symmetric, yet violates both
  identity of indiscernibles and the classical triangle inequality.
  Nevertheless, it obeys the exact *scale-doubled triangle inequality*:
  $
    d_"DOD"^(2 epsilon) (bold(x), bold(z))
      <= d_"DOD"^(epsilon) (bold(x), bold(y)) + d_"DOD"^(epsilon) (bold(y), bold(z))
    wide forall bold(x), bold(y), bold(z) in bb(R)^d.
  $
]

#proof[
  Non-negativity ($d_"DOD"^(epsilon) in {0, dots, d}$), reflexivity
  ($|x_j - x_j| = 0 <= epsilon => d_"DOD"^(epsilon) (bold(x), bold(x)) = 0$),
  and symmetry ($|x_j - y_j| = |y_j - x_j|$) hold directly from the
  definition. Identity of indiscernibles fails whenever
  $0 < norm(bold(x) - bold(y))_infinity <= epsilon$, which yields
  $d_"DOD"^(epsilon) (bold(x), bold(y)) = 0$ with $bold(x) != bold(y)$. To
  refute the classical triangle inequality, take $d = 1$ with $x = 0$,
  $y = 3 epsilon / 4$, and $z = 3 epsilon / 2$; then
  $d_"DOD"^(epsilon) (x, y) = d_"DOD"^(epsilon) (y, z) = 0$, whereas
  $|x - z| = 3 epsilon / 2 > epsilon$ gives $d_"DOD"^(epsilon) (x, z) = 1 > 0 + 0$.
  For the scale-doubled triangle inequality, consider any coordinate
  $j in {1, dots, d}$. Whenever $|x_j - z_j| > 2 epsilon$, the scalar triangle
  inequality $|x_j - y_j| + |y_j - z_j| >= |x_j - z_j| > 2 epsilon$ forces at
  least one of $|x_j - y_j| > epsilon$ or $|y_j - z_j| > epsilon$ to hold.
  Therefore,
  $
    bb(I)(|x_j - z_j| > 2 epsilon)
      <= bb(I)(|x_j - y_j| > epsilon) + bb(I)(|y_j - z_j| > epsilon).
  $
  Summing over $j = 1, dots, d$ proves the exact inequality.
]

#subsec("B", "Escaping High-Dimensional Distance Concentration")

#theorem(title: "Bernstein Contrast Preservation Under Arbitrary Noise Variance", origin: "Original Derivation")[
  Let cluster $C_1$ be localized on an $l$-dimensional subspace $D_1 subset.eq {1, dots, d}$
  such that independent coordinates satisfy:
  + For intra-cluster pairs $(bold(x), bold(y)) in C_1^2$:
    $"Pr"(|x_j - y_j| > epsilon) <= p_"in"$ for $j in D_1$, and
    $"Pr"(|x_j - y_j| <= epsilon) = p_"bg"$ for $j in.not D_1$;
  + For inter-cluster pairs $(bold(x), bold(z))$:
    $"Pr"(|x_j - z_j| <= epsilon) <= p_"bg"$ for all $j in {1, dots, d}$.
  Assume $1 - p_"in" - p_"bg" =: eta > 0$. Then:
  + The expected separation gap satisfies
    $
      Delta_"DOD" := bb(E)[d_"DOD"^(epsilon) (bold(x), bold(z))]
        - bb(E)[d_"DOD"^(epsilon) (bold(x), bold(y))] >= l eta,
    $
    completely independent of the variance or tail heaviness of the $d - l$
    noise coordinates;
  + Setting threshold $theta = d (1 - p_"bg") - l eta / 2$, the pair
    misranking probability obeys the Bernstein tail bound
    $
      "Pr"("misranking") <= 2 exp(- (l^2 eta^2) / (8 V + (4 / 3) l eta)),
    $
    where $V = l p_"in" + d p_"bg"$. In particular, whenever $p_"bg" = O(1 / d)$,
    reliable discrimination holds with $l = Omega(log n)$ regardless of ambient
    dimension $d$.
]

#proof[
  Let $Z_j = bb(I)(|x_j - y_j| > epsilon)$ and $W_j = bb(I)(|x_j - z_j| > epsilon)$.
  By linearity of expectation,
  $bb(E)[d_"DOD"^(epsilon) (bold(x), bold(y))] <= l p_"in" + (d - l)(1 - p_"bg")$,
  while $bb(E)[d_"DOD"^(epsilon) (bold(x), bold(z))] >= d (1 - p_"bg")$.
  Subtracting gives $Delta_"DOD" >= l (1 - p_"bg" - p_"in") = l eta$. Notice
  that heavy-tailed (even infinite-variance Cauchy) noise on $j in.not D_1$
  only decreases $p_"bg"$ and can *only enlarge* $Delta_"DOD"$. For the
  variance of $d_"DOD"^(epsilon) (bold(x), bold(y)) = sum_(j=1)^d Z_j$, each
  independent Bernoulli variable has $"Var"(Z_j) = q_j (1 - q_j)$, bounded by
  $p_"in"$ for $j in D_1$ and by $p_"bg" (1 - p_"bg") <= p_"bg"$ for
  $j in.not D_1$; similarly $"Var"(W_j) <= p_"bg"$ on all $d$ coordinates.
  Thus both sums have variance at most $V = l p_"in" + d p_"bg"$ and centered
  summands bounded in $[-1, 1]$. Applying Bernstein's inequality with
  deviation $t = l eta / 2$ to both $d_"DOD"^(epsilon) (bold(x), bold(y)) >= theta$
  and $d_"DOD"^(epsilon) (bold(x), bold(z)) <= theta$ and summing the two
  events yields
  $
    "Pr"("misranking") <= 2 exp(- (t^2) / (2 V + (2 / 3) t))
      = 2 exp(- (l^2 eta^2) / (8 V + (4 / 3) l eta)).
  $
  When $p_"bg" = O(1 / d)$, $V = l p_"in" + O(1)$, so the exponent scales as
  $-Omega(l)$, driving misranking below $n^(-c)$ as soon as $l = Omega(log n)$
  independently of $d$.
]

#subsec("C", "Exact Chernoff Sample-Size Guarantee for Cluster Discovery")

#theorem(title: "Dimension-Free Chernoff Sample-Size Bound", origin: "Classical — Cited")[
  Suppose $cal(X)$ contains $k$ true clusters $C_1, dots, C_k$, each satisfying
  $|C_r| / n >= alpha_"min" > 0$. Let $S$ be an i.i.d. uniform random sample
  (with replacement) of size $s$. For any required per-cluster sample floor
  $c_"min" >= 1$ and failure probability $delta in (0, 1)$, if
  $
    s >= 1 / alpha_"min" (c_"min" + ln(k / delta)
      + sqrt(ln^2 (k / delta) + 2 c_"min" ln(k / delta))),
  $
  then simultaneously $|S inter C_r| >= c_"min"$ for all $r = 1, dots, k$ with
  probability at least $1 - delta$ @woo2004findit, @guha1998cure.
]

#proof[
  For any cluster $C_r$, let $X_r = |S inter C_r| tilde.op "Bin"(s, alpha_r)$
  with mean $mu_r = s alpha_r >= s alpha_"min" =: mu$. By the multiplicative
  Chernoff lower-tail bound, for $gamma = 1 - c_"min" / mu in (0, 1)$,
  $
    "Pr"(X_r < c_"min") <= "Pr"(X_r <= (1 - gamma) mu_r)
      <= exp(- (mu - c_"min")^2 / (2 mu)).
  $
  Setting $a = c_"min"$ and $b = ln(k / delta)$, requiring
  $(mu - a)^2 / (2 mu) >= b$ is equivalent to the quadratic inequality
  $mu^2 - 2 (a + b) mu + a^2 >= 0$. For $mu > a$, the exact positive root is
  $mu^* = a + b + sqrt(b^2 + 2 a b)$, which satisfies
  $(mu^* - a)^2 = (b + sqrt(b^2 + 2 a b))^2 = 2 b (a + b + sqrt(b^2 + 2 a b)) = 2 b mu^*$
  identically. Hence $s >= mu^* / alpha_"min"$ ensures
  $"Pr"(X_r < c_"min") <= delta / k$ for every $r$. Applying a union bound
  over $r = 1, dots, k$ yields
  $"Pr"(exists r : X_r < c_"min") <= k (delta / k) = delta$.
]

Note that because $sqrt(b^2 + 2 a b) <= a + b$, the cleaner sufficient bound
$s >= (2 / alpha_"min") (c_"min" + ln(k / delta))$ holds unconditionally and
is completely independent of both dataset cardinality $n$ and ambient
dimension $d$.

#subsec("D", "Exact Subspace Recovery via Majority Dimension Voting")

#theorem(title: "Sample Complexity of Exact Subspace Support Recovery", origin: "Original Proof")[
  Let $N_r subset.eq S$ be a candidate cluster neighborhood of size $m = |N_r|$
  around a valid medoid $bold(m)_r in C_r$, containing at most a fraction
  $rho_"cont" in [0, 1/2)$ of contaminant points from outside $C_r$. Suppose
  clean points vote $1$ on each true coordinate $j in D_r$ with probability
  $>= 1 - p_"in"$, while all points vote $1$ on each irrelevant coordinate
  $j in.not D_r$ with probability $<= p_"bg"$. Define
  $mu_"rel" = (1 - rho_"cont")(1 - p_"in")$ and assume $mu_"rel" > p_"bg"$.
  Setting voting threshold $tau_"vote" = (mu_"rel" + p_"bg") / 2$ with margin
  $gamma = (mu_"rel" - p_"bg") / 2 > 0$, the voted subspace $hat(D)_r$
  recovers $D_r$ identically with probability
  $
    "Pr"(hat(D)_r = D_r) >= 1 - d exp(- 2 m gamma^2).
  $
  Hence $m >= (2 gamma^2)^(-1) ln(d / delta')$ neighbors per cluster suffice
  for exact subspace recovery with probability $>= 1 - delta'$.
]

#proof[
  For each coordinate $j in {1, dots, d}$, the empirical vote fraction
  $hat(p)_(r, j) = m^(-1) sum_(bold(x) in N_r) v_j (bold(x), bold(m)_r)$ is an
  average of $m$ independent $[0, 1]$-valued indicators. For any true subspace
  coordinate $j in D_r$, the $(1 - rho_"cont") m$ clean points contribute
  expected vote mass at least $(1 - rho_"cont")(1 - p_"in") = mu_"rel"$, so
  $bb(E)[hat(p)_(r, j)] >= mu_"rel" = tau_"vote" + gamma$. By Hoeffding's
  inequality,
  $"Pr"(j in.not hat(D)_r) = "Pr"(hat(p)_(r, j) < tau_"vote") <= exp(- 2 m gamma^2)$.
  Symmetrically, for any irrelevant coordinate $j in.not D_r$,
  $bb(E)[hat(p)_(r, j)] <= p_"bg" = tau_"vote" - gamma$, so
  $"Pr"(j in hat(D)_r) = "Pr"(hat(p)_(r, j) >= tau_"vote") <= exp(- 2 m gamma^2)$.
  Taking a union bound over all $d$ coordinates gives
  $"Pr"(hat(D)_r != D_r) <= d exp(- 2 m gamma^2)$. Solving
  $d exp(- 2 m gamma^2) <= delta'$ gives $m >= (2 gamma^2)^(-1) ln(d / delta')$,
  scaling only logarithmically in ambient dimension $d$.
]

#subsec("E", "Coordinate Breakdown Robustness vs. Heterogeneous Scale Collapse")

#proposition(title: "Coordinate-Corruption Breakdown and Scale-Sensitivity Limit", origin: "Original Proof")[
  + *Bounded coordinate influence*: Corrupting any $c$ coordinates of a point
    $bold(x)$ by arbitrarily large offsets $Delta_j -> infinity$ alters
    $d_"DOD"^(epsilon) (bold(x), bold(y))$ by at most $c$, preserving exact
    neighborhood membership whenever $c <= d - l_"min" - d_"DOD"^(epsilon) (bold(x), bold(y))$,
    whereas $L_1$ and $L_2$ distances diverge at $c = 1$ (breakdown fraction
    $1 / d$).
  + *Heterogeneous scale failure*: Conversely, under a single global threshold
    $epsilon$, if two unstandardized coordinates have scales
    $sigma_1 << epsilon << sigma_2$, then coordinate $1$ yields
    $p_"bg" -> 1$ (false inclusion in every $hat(D)_r$), while coordinate $2$
    yields $p_"in" -> 1$ (false exclusion from every $hat(D)_r$), driving the
    voting margin $gamma = (mu_"rel" - p_"bg") / 2$ negative.
]

#proof[
  For (1), each indicator $bb(I)(|x_j - y_j| > epsilon) in {0, 1}$ changes by
  at most $1$ regardless of the magnitude of $|x_j - y_j|$. Changing $c$
  coordinates therefore increases $d_"DOD"^(epsilon)$ by at most $c$, whereas
  for any $p >= 1$, $norm(bold(x) - bold(y))_p >= |x_j - y_j| -> infinity$
  upon corrupting even a single coordinate $j$. For (2), if $x_1, y_1$ are
  drawn from any background distribution with variance $sigma_1^2$, Chebyshev's
  inequality gives $"Pr"(|x_1 - y_1| <= epsilon) >= 1 - 2 sigma_1^2 / epsilon^2 -> 1$
  as $sigma_1 / epsilon -> 0$, forcing $p_"bg" -> 1 >= tau_"vote"$ so noise
  axis $1$ is erroneously voted into every subspace. Conversely, on a cluster
  coordinate with $sigma_2 >> epsilon$ and bounded density $f_2 <= C / sigma_2$,
  $"Pr"(|x_2 - m_(r, 2)| <= epsilon) <= 2 C epsilon / sigma_2 -> 0$, forcing
  $1 - p_"in" -> 0 < tau_"vote"$ so signal axis $2$ is discarded.
]

#sec("III", "Complexity, Geometry and Operational Profile")

FINDIT's three-phase architecture decouples quadratic medoid discovery from
dataset size $n$: because sample sizes $|S|$ and $|M|$ depend solely on
$alpha_"min"$, $k$, and $delta$ (Theorem 2), the cluster-forming phase takes
$O(|S| |M| d |cal(E)|)$ time independent of $n$, while the final assignment
scan takes strictly linear $O(n k d)$ time @woo2004findit. On small datasets
where $|S| approx n$, the $O(|S| |M| d)$ pairwise DOD computation introduces
noticeable constant overhead relative to linear grid hashing; on large
datasets ($n >> |S|^2$), FINDIT becomes asymptotically linear in $n$ and $d$
with modest $O(n d + |S| |M|)$ memory footprint.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n k d + |S| |M| d |cal(E)|)$; strictly linear in $n$ once $n >> |S| |M|$],
    [Space Complexity], [$O(n d + |S| |M|)$; single streaming scan suffices for final data assignment],
    [Cluster Geometry], [Disjoint axis-parallel subspace hyperrectangles / $L_1$ balls around medoids],
    [Noise Handling], [Explicit outlier rejection during assignment; DOD caps noisy-axis penalty at $1$],
    [Hyperparameters], [Min cluster fraction $alpha_"min"$, tolerance range $cal(E)$, min subspace dim $l_"min"$, vote $tau_"vote"$],
    [Invariances], [Translation & permutation invariant; uniform scale invariant iff $epsilon$ scales; requires z-scoring],
  ),
  caption: [Operational, geometric, and complexity profile of FINDIT.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *High-Dimensional Biomedical and Clinical Phenotyping*: Electronic health
  record (EHR) cohorts combine $50$--$500$ laboratory biomarkers where distinct
  disease subtypes manifest across different subsets of $5$--$20$ assays while
  the remaining labs exhibit wide patient-to-patient variation. FINDIT's DOD
  metric ignores high-variance irrelevant labs (Theorem 1) while dimension
  voting explicitly outputs the defining clinical biomarker panel $hat(D)_r$
  for each phenotype.
+ *E-Commerce Collaborative User Segmentation*: User-item preference vectors
  in high-dimensional category spaces contain disjoint buyer communities sharing
  tight preferences over distinct product niches. FINDIT discovers both the
  number of niches $k$ and each niche's defining category subset without
  requiring $k$ or $l$ beforehand @woo2004findit.
+ *Industrial Sensor-Fleet Fault Clustering*: Manufacturing lines and aircraft
  turbines monitor hundreds of standardized telemetry channels, where specific
  fault modes perturb only a small subsystem of sensors while unrelated
  channels fluctuate wildly. Proposition 2's coordinate-corruption robustness
  prevents transient glitches on unrelated sensors from distorting fault
  assignments.
+ *High-Dimensional Document and Topic Embedding Partitioning*: Dense term or
  concept feature spaces contain thematic clusters aligned with sparse subsets
  of discriminative topics. Because FINDIT avoids full-space Euclidean norms
  and scales linearly in $n$ during assignment, it handles large corpora
  efficiently.
+ *Fraud Ring and Cyber-Behavior Archetyping*: Financial transaction profiles
  span dozens of velocity, merchant-category, and network indicators where
  different fraud typologies leave signatures in distinct subsets of features;
  FINDIT's explicit outlier rejection separates normal and unclustered noise
  from crisp subspace fraud rings.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

FINDIT degrades in four well-defined regimes: (i) *unstandardized heterogeneous
feature scales* (Proposition 2), which invert the voting margin; (ii) *oblique
or rotated correlated subspaces*, where axis-aligned coordinate indicators
$bb(I)(|x_j - y_j| <= epsilon)$ cannot track diagonal hyperplanes;
(iii) *tiny rare clusters* violating the sampling threshold
$|C_r| / n < alpha_"min"$ (Theorem 2), which may be missed entirely in $M$ and
$S$; and (iv) *overlapping multi-membership clusters*, since FINDIT's
assignment phase enforces a hard disjoint partition.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-D standardized tabular ($d = 20$--$200$)], [Strong: DOD defeats Beyer concentration; discovers $k$ and $D_r$], [Must z-score / IQR-scale all features first],
    [Large-$n$ projected clustering], [Strong: $O(n k d)$ final scan after $n$-independent sampling], [Ensure $|S|$ covers smallest target cluster $alpha_"min"$],
    [Oblique / rotated subspaces], [Poor: axis-aligned DOD misses diagonal correlations], [Use ORCLUS @aggarwal2000orclus or OptiGrid (Alg. 141)],
    [Overlapping subspace patterns], [Weak: hard assignment forces one cluster per point], [Use ENCLUS (Alg. 143) or SUBCLU @kailing2004subclu],
    [Small-$n$ latency-critical tasks], [Moderate: $O(|S| |M| d)$ voting phase dominates at tiny $n$], [Use BANG (Alg. 142) or ENCLUS for low-D tiny $n$],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for FINDIT.],
)

#sec("V", "Structural Relations and Empirical Behavior")

FINDIT belongs to the projected medoid clustering lineage alongside PROCLUS
@aggarwal1999proclus, ORCLUS @aggarwal2000orclus, DOC @procopiuc2002doc,
P3C @moise2006p3c, and COSA @friedman2004cosa, extending $k$-medoids / PAM
@kaufman1990finding to adaptive per-cluster coordinate subspaces while replacing
PROCLUS's fixed-$(k, l)$ local search with DOD neighborhood discovery. The
following proposition formalizes FINDIT's exact algebraic link to PROCLUS and
thresholded Chebyshev / Hamming kernels.

#proposition(title: "Algebraic Reduction Between DOD, Hamming Kernels, and PROCLUS", origin: "Original Proof")[
  Define the coordinate-thresholded similarity kernel
  $K_epsilon (bold(x), bold(y)) = sum_(j=1)^d bb(I)(|x_j - y_j| <= epsilon)$.
  Then:
  + $d_"DOD"^(epsilon) (bold(x), bold(y)) = d - K_epsilon (bold(x), bold(y))$
    equals the exact Hamming distance on pairwise $epsilon$-proximity
    indicators, and the $r$-ball ${bold(y) : d_"DOD"^(epsilon) (bold(x), bold(y)) = 0}$
    coincides with the closed $L_infinity$ (Chebyshev) hypercube
    $B_infinity (bold(x), epsilon)$;
  + Once FINDIT fixes medoids ${bold(m)_r}_(r=1)^k$ and voted subspaces
    ${hat(D)_r}_(r=1)^k$, its assignment phase solves the identical
    per-cluster subspace-normalized $L_1$ objective minimized in the assignment
    step of PROCLUS @aggarwal1999proclus:
    $
      J_"FINDIT" ({C_r}) = sum_(r=1)^k sum_(bold(x) in C_r)
        1 / |hat(D)_r| sum_(j in hat(D)_r) |x_j - m_(r, j)|.
    $
]

#proof[
  For (1), since $bb(I)(|x_j - y_j| > epsilon) = 1 - bb(I)(|x_j - y_j| <= epsilon)$,
  summing over $j = 1, dots, d$ directly gives
  $d_"DOD"^(epsilon) (bold(x), bold(y)) = d - K_epsilon (bold(x), bold(y))$.
  Moreover, $d_"DOD"^(epsilon) (bold(x), bold(y)) = 0$ holds iff
  $|x_j - y_j| <= epsilon$ for every $j$, i.e.,
  $norm(bold(x) - bold(y))_infinity <= epsilon$, which defines
  $B_infinity (bold(x), epsilon)$. For (2), minimizing $J_"FINDIT"$ over
  disjoint assignments $C_1, dots, C_k$ decouples additively across points
  $bold(x) in cal(X)$, so the globally optimal assignment given
  ${(bold(m)_r, hat(D)_r)}$ sends each $bold(x)$ to
  $arg min_r |hat(D)_r|^(-1) sum_(j in hat(D)_r) |x_j - m_(r, j)|$---matching
  PROCLUS's segmental Manhattan assignment rule while replacing PROCLUS's
  user-fixed average dimensionality constraint $sum_r |D_r| = k l$ with
  data-driven dimension voting @kriegel2009subspacesurvey.
]

On the ClusterOverview 4-dataset benchmark (Silhouette @rousseeuw1987silhouettes
and ARI @hubert1985comparing), FINDIT achieves: *Iris* $"Sil" = 0.455$,
$"ARI" = 0.574$; *Customers* $"Sil" = 0.311$; *Multi-Shape* $"Sil" = 0.609$
($"ARI" =$ N/A); *High-D* $"Sil" = 0.264$ ($"ARI" =$ N/A); with a mean runtime
of $2.690$ ms. These empirical measurements reflect every theoretical property
derived above. Most notably, FINDIT achieves the second-highest *High-D*
Silhouette ($0.264$) among Algorithms 141--144---trailing only OptiGrid
($0.287$) and substantially outperforming BANG ($0.168$) and ENCLUS
($0.232$)---confirming Theorem 1's Bernstein contrast preservation under noisy
ambient dimensions. On *Customers* ($"Sil" = 0.311$), FINDIT outperforms both
OptiGrid ($0.232$) and ENCLUS ($0.263$) by isolating correlated spending
subspaces via dimension voting, while its *Iris* score ($"Sil" = 0.455$,
$"ARI" = 0.574$) reflects axis-aligned medoid partitioning. Finally, its
$2.690$ ms mean runtime (higher than the $0.044$--$0.761$ ms grid methods)
corroborates Section III's complexity profile: on small benchmark sizes, the
$O(|S| |M| d)$ pairwise DOD evaluation and threshold-tuning loop dominate,
whereas FINDIT's asymptotic $O(n k d)$ linear scaling pays major dividends at
large $n$.
