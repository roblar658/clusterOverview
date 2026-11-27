#import "../template/ieee.typ": *

#alg_entry(
  id: 134,
  name: "CLIQUE (Clustering in Quest)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "1998",
  paradigm: "Apriori Grid-Based Axis-Parallel Subspace Clustering",
  complexity: $cal(O)(c^k + m k n)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

In high-dimensional data spaces $RR^d$ ($d >> 10$), global pairwise distances concentrate around their mean, eroding nearest-neighbor contrast and causing full-space clustering algorithms to fail @beyer1999nearest. Frequently, however, meaningful clusters exist embedded inside distinct low-dimensional *axis-parallel subspaces* $S subset.eq [d] := {1, dots.h, d}$, while remaining coordinates behave as unstructured noise @kriegel2009subspacesurvey. CLIQUE (CLustering In QUEst) @agrawal1998clique was the first systematic algorithm to combine density-based grid partitioning with bottom-up Apriori search to automatically discover all dense subspaces, their embedded clusters, and concise human-interpretable rule descriptions.

Let $cal(X) = {bold(x)_1, dots.h, bold(x)_n} subset Omega = product_(j=1)^d [a_j, b_j] subset RR^d$ be a finite dataset. CLIQUE fixes two global structural parameters: a per-axis grid resolution $xi in NN_(>= 2)$ and a global density (selectivity) threshold $tau in (0, 1]$. Each coordinate $j in [d]$ is partitioned into $xi$ equal-width half-open intervals $I_(j,1), dots.h, I_(j,xi)$ of width $Delta_j = (b_j - a_j) / xi$ (with the final interval closed on the right).

#definition(title: "Subspace Units, Selectivity, and DNF Cluster Covers")[
  For any non-empty coordinate subset $S = {j_1 < dots.h < j_k} subset.eq [d]$ of cardinality $|S| = k$, a *$k$-dimensional unit* in subspace $S$ is a Cartesian product $u = product_(r=1)^k I_(j_r, t_r)$ with indices $t_r in {1, dots.h, xi}$. Its *selectivity* is the empirical probability mass
  $ sigma(u) = 1 / n |{bold(x) in cal(X) : pi_S (bold(x)) in u}|, $
  where $pi_S : RR^d -> RR^k$ denotes the orthogonal projection onto $S$. Unit $u$ is *dense* if $sigma(u) >= tau$; let $D_S$ denote the set of all dense units in subspace $S$, and $D_k = union.big_(|S|=k) D_S$. Two dense units $u, v in D_S$ are *common-face adjacent* ($u tilde.op_S v$) if they coincide in $k - 1$ coordinates and occupy consecutive intervals ($|t_r - t'_r| = 1$) in the remaining coordinate. A *subspace cluster* $C subset.eq D_S$ is a maximal connected component of the common-face graph $G_S = (D_S, E_S)$. A *minimal DNF description* of $C$ is a minimum-cardinality union of axis-parallel hyper-rectangles (regions) $R_1, dots.h, R_p subset.eq C$ satisfying $union.big_(q=1)^p R_q = C$.
]

CLIQUE executes in three decoupled phases @agrawal1998clique: (1) *Bottom-up subspace discovery* generates candidate $k$-units $C_k$ from dense $(k-1)$-units $D_(k-1)$ via Apriori self-join, scans $cal(X)$ once per level to filter $D_k$, and optionally prunes low-coverage subspaces using a Minimum Description Length (MDL) criterion; (2) *Connected-component extraction* partitions each retained dense subspace $D_S$ into maximal clusters via depth-first search (DFS); and (3) *Minimal description generation* greedily grows maximal hyper-rectangles inside each cluster $C$ and discards redundant regions to produce an irredundant Disjunctive Normal Form (DNF) predicate over attribute ranges.

For MDL subspace pruning at level $k$, let $x_1 >= x_2 >= dots.h >= x_q$ be the sorted coverages $x_(S_j) = sum_(u in D_(S_j)) sigma(u)$ of the $q$ dense $k$-dimensional subspaces. Partitioning the sorted list at split index $i in {1, dots.h, q - 1}$ into selected set $I = {1, dots.h, i}$ and pruned set $P = {i + 1, dots.h, q}$ with means $mu_I (i)$ and $mu_P (i)$ incurs total encoding cost
$ L_"MDL" (i) &= log_2 mu_I (i) + sum_(j=1)^i log_2 |x_j - mu_I (i)| \
  &quad + log_2 mu_P (i) + sum_(j=i+1)^q log_2 |x_j - mu_P (i)|, $
where CLIQUE retains only subspaces in $I$ at $i^* = arg min_i L_"MDL" (i)$.

#algbox(
  title: "CLIQUE Bottom-Up Subspace Clustering",
  inputs: [Dataset $cal(X) subset RR^d$, grid intervals $xi$, selectivity threshold $tau in (0, 1]$, optional MDL pruning flag.],
  outputs: [Collection of subspace clusters ${(S_m, C_m, cal(R)_m)}$ with DNF hyper-rectangle covers $cal(R)_m$.],
)[
  + *Level-1 Scan:* Partition each axis $j in [d]$ into $xi$ bins; scan $cal(X)$ to find 1D dense units $D_1 = {u : sigma(u) >= tau}$.
  + *Apriori Lattice Search:* Set $k <- 2$. While $D_(k-1) != emptyset$:
    - *Candidate Join:* Form $C_k$ by joining pairs $u, v in D_(k-1)$ sharing their first $k - 2$ dimensions and intervals, with $j_(k-1)(u) < j_(k-1)(v)$.
    - *Subset Pruning:* Discard any $w in C_k$ having a $(k-1)$-face projection $w_(-q) in.not D_(k-1)$.
    - *Selectivity Scan:* Scan $cal(X)$ once to retain $D_k = {w in C_k : sigma(w) >= tau}$.
    - *MDL Pruning (optional):* Sort subspaces by coverage $x_S$; retain $I = {1, dots.h, i^*}$ minimizing $L_"MDL" (i)$.
    - Increment $k <- k + 1$.
  + *Cluster Extraction:* For each retained subspace $S$ with $D_S != emptyset$, run DFS on $G_S = (D_S, E_S)$ to identify maximal connected components $C$.
  + *Minimal DNF Covering:* For each cluster $C$, greedily expand uncovered dense units into maximal axis-aligned hyper-rectangles $R subset.eq C$, then iteratively remove any region $R$ fully covered by $union.big_(R' != R) R'$.
  + *Return* all subspace clusters and their irredundant DNF covers.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five rigorous theorems analyzing CLIQUE's downward-closure completeness, the dimensionality-threshold squeeze caused by a global $tau$, exact scale invariance and recall loss under MDL pruning, exponential lattice lower bounds plus NP-hardness of minimal DNF covers, and geometric invariance properties.

#subsec("A", "Downward Closure and Apriori Completeness")

#theorem(title: "Monotonicity of Density and Soundness of Apriori Pruning", origin: "Classical — Cited")[
  If a $k$-dimensional unit $u = product_(r=1)^k I_(j_r, t_r)$ in subspace $S$ is dense ($sigma(u) >= tau$), then every $(k-1)$-dimensional projection $u_(-q) = product_(r != q) I_(j_r, t_r)$ onto $S without {j_q}$ ($q in {1, dots.h, k}$) is also dense @agrawal1998clique. Consequently, without MDL pruning, CLIQUE's candidate generation and subset-pruning step is sound and complete: $D_k$ equals the exact set of all dense $k$-dimensional units in $Omega$.
]

#proof[
  For any observation $bold(x) in cal(X)$ and any coordinate index $q in {1, dots.h, k}$, membership $pi_S (bold(x)) in u$ requires $x_(j_r) in I_(j_r, t_r)$ for all $r = 1, dots.h, k$, which strictly implies $x_(j_r) in I_(j_r, t_r)$ for all $r != q$, i.e., $pi_(S without {j_q}) (bold(x)) in u_(-q)$. Hence
  $ {bold(x) in cal(X) : pi_S (bold(x)) in u} subset.eq {bold(x) in cal(X) : pi_(S without {j_q}) (bold(x)) in u_(-q)}. $
  Dividing cardinalities by $n$ gives $sigma(u_(-q)) >= sigma(u) >= tau$, establishing downward closure.

  We prove completeness of $D_k$ by induction on $k >= 1$. For $k = 1$, exhaustive binning finds all dense 1-units. Assume $D_(k-1)$ is exact and let $w = product_(r=1)^k I_(j_r, t_r)$ be any true dense $k$-unit with $j_1 < dots.h < j_k$. By monotonicity, both $u = w_(-k)$ (dimensions $j_1, dots.h, j_(k-1)$) and $v = w_(-(k-1))$ (dimensions $j_1, dots.h, j_(k-2), j_k$) belong to $D_(k-1)$. Because $u$ and $v$ coincide on their first $k - 2$ dimensions and satisfy $j_(k-1) < j_k$, the self-join explicitly generates $w in C_k$. Moreover, every other $(k-1)$-face $w_(-q)$ belongs to $D_(k-1)$ by monotonicity, so $w$ survives subset pruning and passes the exact selectivity scan $sigma(w) >= tau$. Soundness is immediate because every unit in $D_k$ is verified directly against $cal(X)$.
]

#subsec("B", "Coverage Monotonicity and the Global-Threshold Dimensionality Squeeze")

A fundamental mathematical limitation of CLIQUE is enforcing one fixed global threshold $tau$ across all dimensionalities $k = 1, dots.h, d$. We formalize the exact trade-off between detecting high-dimensional clusters and rejecting low-dimensional background noise.

#theorem(title: "Coverage Monotonicity and Sensitivity-Specificity Squeeze", origin: "Original Derivation")[
  Let $C(S) = sum_(u in D_S) sigma(u)$ denote the total empirical coverage of subspace $S$.
  + *Chain Monotonicity:* For any nested subspaces $S subset S' subset.eq [d]$, empirical coverage is non-increasing: $C(S') <= C(S)$.
  + *Dimensionality Squeeze:* Suppose $cal(X)$ consists of a cluster of mass fraction $alpha in (0, 1)$ uniformly spread across $m_0 >= 2$ grid intervals along each of $k^*$ relevant dimensions $S^*$ (so $m_0^(k^*)$ cluster units), plus uniform background noise of mass $1 - alpha$ across all $xi$ intervals on every axis. Ignoring finite-sample fluctuations, simultaneously detecting the $k^*$-dimensional cluster in $S^*$ and rejecting spurious 1D noise units on irrelevant axes requires
    $ (1 - alpha) / xi < tau <= alpha / m_0^(k^*), $
    which is *mathematically impossible* whenever $k^* > log_(m_0) (alpha xi / (1 - alpha))$.
]

#proof[
  For (1), let $cal(X)_D (S) = {bold(x) in cal(X) : exists u in D_S, pi_S (bold(x)) in u}$. Because units in $D_S$ are pairwise disjoint, $C(S) = |cal(X)_D (S)| / n$. Take any $bold(x) in cal(X)_D (S')$, so $pi_(S') (bold(x)) in w$ for some dense unit $w in D_(S')$. By repeated application of Theorem II.A along the inclusion $S subset S'$, the projection $u = pi_S (w)$ satisfies $sigma(u) >= sigma(w) >= tau$, hence $u in D_S$ and $pi_S (bold(x)) in u$. Thus $cal(X)_D (S') subset.eq cal(X)_D (S)$, proving $C(S') <= C(S)$.

  For (2), inside the relevant subspace $S^*$ of dimension $k^*$, the $alpha n$ cluster points distribute evenly across $m_0^(k^*)$ congruent units, while background noise contributes $(1 - alpha) / xi^(k^*)$ per unit. Assuming negligible background contribution inside a $k^*$-unit ($xi >> m_0$), each cluster unit achieves selectivity $sigma(u) approx alpha / m_0^(k^*)$, so detecting the cluster requires $tau <= alpha / m_0^(k^*)$. Conversely, on any irrelevant coordinate $j in.not S^*$, background noise alone assigns selectivity $(1 - alpha) / xi$ to every 1D interval (plus $alpha / xi$ from the cluster if uniformly projected). Preventing all $xi$ bins on every irrelevant axis from turning dense requires $tau > (1 - alpha) / xi$. Combining both inequalities gives $(1 - alpha) / xi < alpha / m_0^(k^*)$, equivalently $m_0^(k^*) < alpha xi / (1 - alpha)$, or $k^* < log_(m_0) (alpha xi / (1 - alpha))$. Above this critical dimension, *no single value of $tau$* can recover the $k^*$-dimensional cluster without flooding $D_1$ with dense noise bins across all $d$ coordinates.
]

#subsec("C", "Exact Scale Invariance and Recall Loss of MDL Subspace Pruning")

#proposition(title: "Properties and Incompleteness of MDL Pruning", origin: "Original Proof")[
  Let $x_1 >= dots.h >= x_q > 0$ be distinct subspace coverages at level $k$ with optimal MDL cut $i^* = arg min_(1 <= i < q) L_"MDL" (i)$.
  + *Exact Rescaling Invariance:* For every positive constant $gamma > 0$, rescaling all coverages $x_j |-> gamma x_j$ shifts $L_"MDL" (i)$ by an $i$-independent constant $(q + 2) log_2 gamma$, leaving the optimal split $i^*$ strictly unchanged.
  + *Loss of Completeness:* There exists a valid dataset with two disjoint true clusters in subspaces $S_1, S_2$ ($|S_1| = |S_2| = 2$) where MDL pruning discards $S_2$ even though $S_2$ contains dense units well above $tau$, breaking Theorem II.A's zero-false-negative guarantee.
]

#proof[
  For (1), linearity of the arithmetic mean gives $mu'_I (i) = gamma mu_I (i)$ and $mu'_P (i) = gamma mu_P (i)$. Substituting into $L_"MDL" (i)$ yields
  $ L'_"MDL" (i) &= log_2 (gamma mu_I (i)) + sum_(j=1)^i log_2 (gamma |x_j - mu_I (i)|) \
    &quad + log_2 (gamma mu_P (i)) + sum_(j=i+1)^q log_2 (gamma |x_j - mu_P (i)|) \
    &= L_"MDL" (i) + (2 + i + (q - i)) log_2 gamma \
    &= L_"MDL" (i) + (q + 2) log_2 gamma. $
  Because $(q + 2) log_2 gamma$ is constant with respect to $i$, $arg min_i L'_"MDL" (i) = arg min_i L_"MDL" (i) = i^*$.

  For (2), suppose a dataset has one large cluster embedded in a 4D subspace generating $binom(4, 2) = 6$ level-2 subspaces with tightly clustered coverages $x_1, dots.h, x_6 in [0.700, 0.705]$ (so $|x_j - mu_I (6)| <= 0.0025$), one smaller true cluster of coverage $x_7 = 0.12 >= tau = 0.05$ in a disjoint 2D subspace $S_2$, and two low-coverage subspaces $x_8 = 0.08$, $x_9 = 0.06$. At split $i = 6$, every residual in $I$ has negative bit-length contribution $log_2 (0.0025) approx -8.64$, whereas moving $S_2$ into $I$ at $i = 7$ inflates all seven residuals to $|x_j - mu_I (7)| >= 0.08$, raising $L_"MDL" (7)$ by more than $35$ bits over $L_"MDL" (6)$. Thus $i^* = 6$ prunes $S_2$, causing zero recall on the second cluster and all its super-spaces.
]

#subsec("D", "Exponential Lattice Lower Bound and DNF Cover Hardness")

#theorem(title: "Unconditional Exponential Subspace Lower Bound and Cover Complexity", origin: "Original Proof")[
  + *Exponential Lattice Lower Bound:* If $cal(X)$ contains at least one dense unit $u^*$ in a $k$-dimensional subspace $S^*$, then unpruned CLIQUE must generate, store, and verify at least $2^k - 1$ distinct non-empty dense subspaces, imposing an unconditional $Omega(2^k)$ lower bound on time and output cardinality regardless of $n$.
  + *Hardness of Optimal DNF Covering:* Computing a minimum-cardinality union of axis-parallel hyper-rectangles covering an arbitrary connected union of dense grid units $C subset.eq {1, dots.h, xi}^k$ is NP-hard for every $k >= 2$, whereas CLIQUE's greedy maximal-region cover with redundant removal terminates in polynomial time in $|C|$ and always outputs an inclusion-minimal (irredundant) cover.
]

#proof[
  For (1), let $u^* in D_(S^*)$ with $|S^*| = k$. By Theorem II.A, for every non-empty subset $T subset.eq S^*$, the projection $pi_T (u^*)$ satisfies $sigma(pi_T (u^*)) >= sigma(u^*) >= tau$, so $D_T != emptyset$. Because $S^*$ has exactly $sum_(j=1)^k binom(k, j) = 2^k - 1$ non-empty subsets and CLIQUE enumerates all dense units level by level from $j = 1$ to $k$, the algorithm must create at least $2^k - 1$ dense units across $2^k - 1$ distinct subspaces before discovering $u^*$, yielding an $Omega(2^k)$ lower bound (matching the classical $cal(O)(c^k + m k n)$ upper bound @agrawal1998clique).

  For (2), for $k = 2$, any orthogonal polygon with holes on an integer grid is a connected union of unit squares $C subset.eq {1, dots.h, xi}^2$, and finding a minimum number of axis-aligned rectangles covering an orthogonal polygon with holes is classically NP-hard (by reduction from Planar 3-SAT); embedding the 2D grid into $k > 2$ dimensions by fixing coordinates $3, dots.h, k$ to a single interval preserves NP-hardness for all $k >= 2$. In CLIQUE's heuristic, each greedy step starts from an uncovered unit $u in C$ and expands maximally along each of the $k$ axes in $cal(O)(k |C|)$ time; since each step covers at least one new unit, at most $|C|$ maximal regions are created, taking $cal(O)(k |C|^2)$ time. The subsequent removal pass checks whether each region $R$ is covered by $union.big_(R' != R) R'$ in $cal(O)(|C|^2)$ time and deletes it if so; upon termination, no remaining region can be deleted without exposing at least one unit of $C$, guaranteeing irredundancy.
]

#subsec("E", "Axis-Affine Invariance vs. Oblique Rotation Fragility")

#proposition(title: "Axis-wise Affine Invariance and Oblique Attenuation", origin: "Original Proof")[
  CLIQUE's dense units, connected clusters, and combinatorial DNF descriptions are strictly invariant under data-order permutations $pi in S_n$ and arbitrary strictly increasing per-axis affine transformations $x_(i,j) |-> alpha_j x_(i,j) + beta_j$ ($alpha_j > 0$). However, rotating a 1D axis-aligned line cluster in $RR^2$ by $pi / 4$ reduces its peak 2D unit selectivity by factor $Theta(xi)$.
]

#proof[
  Under $tilde(x)_(i,j) = alpha_j x_(i,j) + beta_j$ ($alpha_j > 0$), domain endpoints transform identically as $tilde(a)_j = alpha_j a_j + beta_j$ and $tilde(b)_j = alpha_j b_j + beta_j$, so bin width scales to $tilde(Delta)_j = alpha_j Delta_j$. Thus the grid index
  $ floor((tilde(x)_(i,j) - tilde(a)_j) / tilde(Delta)_j) = floor((alpha_j (x_(i,j) - a_j)) / (alpha_j Delta_j)) = floor((x_(i,j) - a_j) / Delta_j) $
  is preserved identically for every point and coordinate. Since selectivities $sigma(u)$ are symmetric counting functions, both permutation and per-axis positive affine invariance follow immediately.

  Conversely, place $n$ points uniformly along the vertical segment ${0.5} times [0, 1] subset [0, 1]^2$. In 2D, every point falls into one vertical column of $xi$ units, each with selectivity $1 / xi$, while its 1D $x_1$-projection concentrates mass $1$ into a single bin. Rotating by $pi / 4$ maps the segment onto the diagonal $x_1 = x_2$, intersecting $xi$ diagonal 2D units (selectivity $1 / xi$ each) while spreading *both* 1D marginal projections uniformly across $xi$ bins of selectivity $1 / xi$ — reducing peak 1D selectivity by exact factor $xi$ and causing the entire subspace to vanish for any $tau in (1 / xi, 1]$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

For fixed maximal dense subspace dimension $k$, CLIQUE makes exactly $k$ sequential passes over the dataset $cal(X)$ — one pass per Apriori level — requiring $cal(O)(m k n)$ data-scan operations where $m$ is the maximum number of candidate subspaces evaluated at a level, plus $cal(O)(c^k)$ ($c > 1$) combinatorial candidate-generation overhead @agrawal1998clique. Consequently, runtime scales *linearly* with sample size $n$ (making CLIQUE highly scalable for large $n$ in low-$k$ regimes) yet *exponentially* with the intrinsic dimension $k$ of the highest-dimensional dense cluster (Theorem II.D). Memory scales as $cal(O)(|C_k|)$ to store hash tables of active dense units at level $k$.

Geometrically, CLIQUE produces *overlapping subspace clusterings*: a single point $bold(x)_i$ may simultaneously belong to distinct clusters in multiple subspaces $S_1, S_2$ (as well as every lower-dimensional projection of a dense $k$-cluster). Within each subspace $S$, clusters are arbitrary connected unions of congruent hypercubes, capturing non-convex manifolds while emitting explicit human-readable DNF range rules.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$cal(O)(c^k + m k n)$: linear in $n$, exponential in max cluster dimension $k$],
    [Space Complexity], [$cal(O)(|C_k| + |D_k|)$ active hash-table units across adjacent levels],
    [Database Passes], [Exactly $k$ sequential scans over $cal(X)$ (one per Apriori level)],
    [Cluster Overlap], [Overlapping across subspaces; disjoint connected components within each $S$],
    [Invariances], [Permutation-invariant and per-axis positive affine-invariant; not rotation-invariant],
    [Interpretability], [High: outputs minimal DNF unions of axis-parallel interval rules],
    [Primary Hyperparams], [Grid resolution $xi$ (bins/axis), global selectivity threshold $tau$],
    [Key Limitations], [$Omega(2^k)$ blowup for large $k$; fixed-$tau$ dimensionality bias; grid boundaries],
  ),
  caption: [Computational, geometric, and operational profile of CLIQUE @agrawal1998clique.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

*1. Gene Expression and Multi-Omics Subspace Discovery:* In transcriptomics and proteomics, biological pathways co-regulate subsets of genes only under specific subsets of experimental conditions or tissue types. CLIQUE's overlapping subspace paradigm naturally discovers multiple co-active regulatory modules where the same gene participates in distinct pathways across different feature subsets, while its DNF intervals translate directly into thresholded expression rules @kriegel2009subspacesurvey.

*2. High-Dimensional Customer Segmentation and Targeted Marketing:* Enterprise CRM tables combine hundreds of demographic, behavioral, and transactional features, yet actionable consumer cohorts typically cohere along 3–5 specific attributes (e.g., age, tenure, and two product categories). CLIQUE's per-axis affine invariance (Proposition II.E) eliminates ad-hoc feature scaling between currencies and ages, and its minimal DNF descriptions compile straight into SQL marketing filters @agrawal1998clique.

*3. Network Intrusion Detection and Security Forensics:* Packet-flow telemetry records dozens of protocol, byte-count, and timing descriptors where distinct attack signatures localize within small axis-parallel feature subsets (e.g., SYN rate plus packet-size variance). Bottom-up Apriori search exposes multi-feature attack modes while ignoring irrelevant traffic fields, and the resulting hyper-rectangle rules serve directly as interpretable firewall signatures.

*4. Financial Factor and Market-Regime Mining:* Quantitative analysts track dozens of macroeconomic and technical indicators where asset co-movements emerge transiently inside low-dimensional factor subspaces. CLIQUE detects overlapping factor regimes without forcing every asset into a single global partition, accommodating multi-factor exposure naturally.

*5. Cheminformatics and QSAR Descriptor Screening:* Molecular libraries represented by hundreds of physicochemical descriptors often exhibit bioactivity cliffs confined to narrow ranges of 2–4 key descriptors (e.g., lipophilicity, polar surface area, molecular weight). CLIQUE's grid density search identifies high-activity descriptor pockets and outputs compact bounding-box rules useful for medicinal chemistry optimization.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CLIQUE exhibits four canonical failure modes: (i) *exponential blowup* when intrinsic cluster dimensionality $k > 10$ (Theorem II.D); (ii) *dimensionality-threshold collapse* when global $tau$ either misses high-$k$ clusters or floods $D_1$ with noise (Theorem II.B); (iii) *oblique / correlated subspace blindness* when clusters lie on rotated hyperplanes (Proposition II.E); and (iv) *redundant output proliferation*, emitting up to $2^k - 1$ lower-dimensional shadow clusters for every $k$-dimensional cluster.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [CRM / Interpretability Rules], [Exact DNF interval rules; per-axis affine invariance], [Shadow duplicates; filter maximal subspaces via SCHISM / INSCY @assent2008inscy],
    [Omics / Pathway Subspaces], [Overlapping subspace membership matches multi-role genes], [Large $k > 10$ stalls; use PROCLUS @aggarwal1999proclus or FIRES @kriegel2005fires],
    [Moderate-D Tabular ($d <= 20$)], [Fast linear-$n$ scans when true $k <= 5$], [Fixed $xi$ splits clusters; use adaptive MAFIA @goil1999mafia],
    [Correlated / Rotated Manifolds], [Poor: axis-aligned bins dilute diagonal density by $xi$], [Use ORCLUS @aggarwal2000orclus or SSC],
    [High-D Dense Clusters ($k >= 12$)], [Poor: $Omega(2^k)$ lattice explosion + $tau$ squeeze], [Use PROCLUS @aggarwal1999proclus or SUBCLU @kailing2004subclu],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for CLIQUE.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CLIQUE serves as the canonical ancestor of the entire bottom-up axis-parallel subspace clustering taxonomy @kriegel2009subspacesurvey @aggarwal2013data.

#proposition(title: "Taxonomy Reduction Across Subspace Clustering Families", origin: "Original Proof")[
  Let $cal(H)_"CLIQUE" (xi, tau)$ denote CLIQUE's Apriori subspace lattice. Three landmark successors arise as exact structural generalizations of $cal(H)_"CLIQUE"$:
  + *Adaptive Grid Generalization (MAFIA @goil1999mafia):* Replacing equal-width intervals $Delta_j$ with variable-width histogram-merged intervals $I_(j,t)$ preserves Theorem II.A's monotonicity while tightening the denominator $m_0^(k^*)$ in Theorem II.B toward $1$.
  + *Entropy-Criterion Generalization (ENCLUS @cheng1999enclus):* Replacing density selectivity $sigma(u) >= tau$ with upper-bounded Shannon grid entropy $H(S) = -sum_u sigma(u) log_2 sigma(u) <= omega$ satisfies *upward closure* ($S subset S' => H(S') >= H(S)$), enabling identical Apriori pruning of high-entropy (uniform) subspaces while additionally prioritizing high-interest-gain subspaces.
  + *Grid-Free Density Limit (SUBCLU @kailing2004subclu):* Letting grid boundaries slide continuously to $epsilon$-balls around every point replaces grid connectivity with DBSCAN density connectivity @ester1996dbscan, which retains downward closure ($N_epsilon^(S') (x) subset.eq N_epsilon^S (x)$) and eliminates grid-phase splitting at the cost of $cal(O)(n^2)$ per subspace.
]

#proof[
  For (1), any disjoint partition of axis $j$ into arbitrary intervals ${I_(j,t)}$ satisfies $pi_S (bold(x)) in product I_(j_r, t_r) => pi_(S without {j_q}) (bold(x)) in product_(r != q) I_(j_r, t_r)$, preserving exact monotonicity; merging adjacent dense histogram bins collapses a cluster spanning $m_0$ uniform bins into a single adaptive interval ($m_0 -> 1$), eliminating exponential dilution $m_0^(-k^*)$.

  For (2), because $(pi_S, pi_(S' without S))$ refines the marginal grid distribution on $S$, the chain rule of discrete Shannon entropy gives $H(S') = H(S) + H(S' without S | S) >= H(S)$ (since conditional discrete entropy is non-negative). Hence $H(S') <= omega => H(S) <= omega$, proving downward closure of low-entropy subspaces for Apriori search.

  For (3), for any $S subset S'$, $norm(pi_(S') (x - y))_2 >= norm(pi_S (x - y))_2$, so every $epsilon$-neighborhood in $S'$ is contained in the $epsilon$-neighborhood in $S$: $N_epsilon^(S') (x) subset.eq N_epsilon^S (x)$. Therefore every DBSCAN core point and density-connected cluster in $S'$ remains a core point and density-connected in $S$, transferring Apriori completeness verbatim to grid-free density clustering.
]

Contrast this bottom-up overlapping lattice with top-down projected clustering (PROCLUS @aggarwal1999proclus, Alg. 136, and DOC @procopiuc2002doc), which sacrifices exhaustive subspace enumeration to obtain a disjoint partition in polynomial time in $d$, and with statistical cell-based methods such as P3C @moise2006p3c that replace the global $tau$ with Poisson $p$-values to counteract Theorem II.B's dimensionality bias @hastie2009elements @jain2010data.

*Empirical Benchmark Behavior:* On the ClusterOverview 4-dataset benchmark evaluated via Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing, CLIQUE obtains *Iris* $"Sil" = 0.611$, $"ARI" = 0.543$; *Customers* $"Sil" = 0.176$; *Multi-Shape* $"Sil" = 0.415$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = "N/A"$ ($"ARI" = "N/A"$), with a mean runtime of *$0.145$ ms*. Every single metric directly corroborates Section II's theoretical analysis:
- *Low-D Subspace Peak on Iris:* CLIQUE achieves the *highest* Silhouette ($0.611$) and *highest* ARI ($0.543$) among all four algorithms (Algs. 133–136). In Iris ($d = 4$), species separation concentrates strongly in the 2D petal-length / petal-width subspace; CLIQUE's bottom-up subspace search isolates this exact axis-aligned subspace without distortion from noisier sepal dimensions.
- *Continuous-Density Weakness on Customers & Multi-Shape:* On *Customers* ($"Sil" = 0.176$, lowest of the quartet) and *Multi-Shape* ($"Sil" = 0.415$), coarse fixed-$xi$ binning and a single global $tau$ struggle with smooth skewed spending gradients and curved boundaries.
- *High-D Collapse Confirming Theorem II.B:* On *High-D*, Silhouette is `"N/A"` — a textbook empirical manifestation of Theorem II.B's sensitivity-specificity squeeze, where fixed-$tau$ selectivity dilution ($alpha / m_0^k$) produces either zero valid multi-cluster partitions or a single degenerate cluster, underscoring why polynomial projected methods like PROCLUS ($"Sil" = 0.294$) dominate genuinely high-dimensional benchmarks.
