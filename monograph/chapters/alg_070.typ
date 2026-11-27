#import "../template/ieee.typ": *

#alg_entry(
  id: 70,
  name: "HCA-DBSCAN: Hierarchical Linkage on Density Cores",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2016",
  paradigm: "Two-Stage Density-Core Extraction and Agglomerative Linkage",
  complexity: $O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Density-based clustering via DBSCAN @ester1996dbscan @sander1998gdbscan
excels at isolating arbitrary non-convex shapes and rejecting low-density
noise, yet suffers from three structural limitations: (i) a single
global neighborhood scale $epsilon$ cannot resolve clusters of markedly
different densities without either fracturing sparse clusters or merging
adjacent dense modes @ankerst1999optics @campello2013hdbscan; (ii)
classical DBSCAN provides a flat partition rather than a multi-scale
hierarchy; and (iii) border points adjacent to multiple cores exhibit
input-order ambiguity @schubert2017dbscanrevisited. Conversely, classical
Hierarchical Agglomerative Clustering (*HCA*) @ward1963hierarchical
@mullner2011modern yields rich multi-scale dendrograms but incurs
$O(n^2)$ linkage costs and degrades severely under noisy bridges.

*HCA-DBSCAN* unifies both paradigms in a principled three-stage hybrid
architecture: (1) extract fine-grained *density cores* at a conservative
neighborhood radius $epsilon$ to purge noise and compress high-density
regions; (2) construct a full agglomerative dendrogram over the $m << n$
extracted density cores using a Lance--Williams linkage recurrence
@lance1967general; and (3) assign border points deterministically to
their nearest core's macro-cluster while retaining isolated noise.

Let $(X, d)$ be a finite metric space with $X = {x_1, dots, x_n}$,
neighborhood radius $epsilon > 0$, and minimum density integer
$mu in bb(N)$ ($"MinPts"$). Define the closed $epsilon$-ball
$N_epsilon (x) = {y in X : d (x, y) <= epsilon}$ and partition $X$ into:
$
  X_"core" &= {x in X : abs(N_epsilon (x)) >= mu}, \
  X_"border" &= {x in X without X_"core" :
    N_epsilon (x) inter X_"core" != emptyset}, \
  X_"noise" &= X without (X_"core" union X_"border").
$

#definition(title: "Density-Core Graph and HCA-DBSCAN Hierarchy")[
  Let $G_epsilon = (X_"core", E_epsilon)$ be the *core-connectivity
  graph* with edge $(x, y) in E_epsilon$ iff $d (x, y) <= epsilon$.
  The *density cores* $cal(M) = {M_1, dots, M_m}$ are the connected
  components of $G_epsilon$. Given a base inter-core dissimilarity
  $D_0 : cal(M) times cal(M) -> bb(R)_(>= 0)$ and Lance--Williams
  parameters $(alpha_A, alpha_B, beta, gamma)$, HCA-DBSCAN builds a
  dendrogram $cal(T)_(cal(M))$ over $cal(M)$ by iteratively merging
  $(A^*, B^*) = arg min_(A != B) D (A, B)$ and updating
  $
    D (A union B, C) &= alpha_A D (A, C) + alpha_B D (B, C) \
      &quad + beta D (A, B) + gamma |D (A, C) - D (B, C)|.
  $
  Cutting $cal(T)_(cal(M))$ into $K$ macro-clusters
  $cal(C) = {C_1, dots, C_K}$ induces a partition of $X_"core"$,
  extended to $X_"border"$ via deterministic nearest-core assignment
  $c (y) = c (arg min_(x in X_"core") d (y, x))$.
]

Operating at a tight radius $epsilon$ ensures that distinct high-density
modes remain cleanly separated into distinct components $M_a$, while any
over-fragmentation of lower-density regions is gracefully reassembled in
Stage 2 along the core dendrogram—achieving multi-density sensitivity
with noise immunity and quadratic compression from $n$ points to $m$
cores.

#algbox(
  title: "HCA-DBSCAN (Hierarchical Linkage on Density Cores)",
  inputs: [Dataset $(X, d)$; core parameters $(epsilon, mu)$; linkage
    rule $(alpha, beta, gamma)$; cut criterion $K$ or gap threshold],
  outputs: [Core dendrogram $cal(T)_(cal(M))$, cluster labels on
    $X_"core" union X_"border"$, and noise set $X_"noise"$],
)[
  + *Stage 1 (Core Extraction):* For each $x in X$, query
    $N_epsilon (x)$ (via spatial index in low $d$, or distance scan).
    Identify $X_"core" = {x : abs(N_epsilon (x)) >= mu}$, $X_"border"$,
    and $X_"noise"$.
  + *Core Components:* Extract connected components
    $cal(M) = {M_1, dots, M_m}$ of $G_epsilon = (X_"core", E_epsilon)$
    via BFS/DFS in $O(abs(E_epsilon))$.
  + *Stage 2 (Core Linkage):* Compute pairwise inter-core
    dissimilarities $D (M_a, M_b)$ (Single, Complete, Average, or
    weighted Ward over core centroids $(abs(M_a), bar(x)_(M_a))$). Run
    agglomerative clustering (via NN-chain for reducible linkages in
    $O(m^2)$) to form $cal(T)_(cal(M))$.
  + *Stage 3 (Cut & Deterministic Border Assignment):* Cut
    $cal(T)_(cal(M))$ into $K$ macro-clusters ${C_1, dots, C_K}$. For
    each $y in X_"border"$, assign $y$ to the cluster of
    $x^* (y) = arg min_(x in X_"core") d (y, x)$; label $X_"noise"$ as
    noise ($-1$).
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove four core mathematical properties of HCA-DBSCAN: an exact
equivalence to $epsilon$-truncated Single Linkage on cores paired with
unbounded noise-bridge immunity, necessary and sufficient reducibility
conditions guaranteeing ultrametricity and strict permutation invariance,
classical density level-set consistency with multi-density separation,
and an exact Huygens--Steiner variance equivalence for core-compressed
Ward linkage.

#subsec("A", "Truncated Single-Linkage Equivalence and Bridge Immunity")

#theorem(
  title: "Core-Linkage Equivalence and Unbounded Noise-Bridge Immunity",
  origin: "Original Proof",
)[
  Let $cal(M) = {M_1, dots, M_m}$ be the density cores of $(X, d)$ at
  $(epsilon, mu)$.
  + When Stage 2 uses Single Linkage $D_"SL" (M_a, M_b) = min_(x in M_a,
    y in M_b) d (x, y)$, the resulting core dendrogram $cal(T)_(cal(M))$
    is isomorphic to applying Single Linkage directly to $X_"core"$ and
    contracting all merges of height $<= epsilon$.
  + *Noise-Bridge Immunity:* Suppose two core sets $A, B subset.eq
    X_"core"$ are separated by $Delta = min_(x in A, y in B) d (x, y) >
    2 epsilon$, and linked by a sparse noise bridge $Z = {z_1, dots,
    z_L}$ with consecutive spacing $d (z_r, z_(r + 1)) <= delta << Delta$
    such that $abs(N_epsilon (z_r)) < mu$ for every $z_r in Z$. Then
    standard Single Linkage merges $A$ and $B$ at height $<= delta$,
    whereas HCA-DBSCAN preserves $A$ and $B$ as separate clusters up to
    merge height $Delta$—an unbounded robustness gain $Delta / delta ->
    infinity$.
]

#proof[
  (1) Consider Kruskal's MST construction on the complete graph over
  $X_"core"$. Every edge $(x, y)$ inside $G_epsilon$ has weight
  $<= epsilon$, whereas every edge joining distinct connected components
  $M_a != M_b$ of $G_epsilon$ has weight strictly greater than $epsilon$
  (otherwise $M_a$ and $M_b$ would belong to the same connected
  component). Sorting all pairwise edges in non-decreasing weight order
  therefore processes *all* edges of weight $<= epsilon$ strictly before
  any edge of weight $> epsilon$. After processing all edges of weight
  $<= epsilon$, the connected components of the Kruskal forest are
  precisely $cal(M) = {M_1, dots, M_m}$. Subsequent Kruskal merges over
  edges $> epsilon$ select at each step the lightest edge connecting two
  active unions of components, which is identically Single Linkage on
  $(cal(M), D_"SL")$.

  (2) On $X = A union B union Z$, the path $a - z_1 - dots.c - z_L - b$
  (with $d (a, z_1) <= delta$ and $d (z_L, b) <= delta$) has bottleneck
  edge weight $<= delta$, so standard Single Linkage joins $A$ and $B$
  at merge height $<= delta$. Under HCA-DBSCAN, because
  $abs(N_epsilon (z_r)) < mu$ for every $r$, we have $Z inter X_"core"
  = emptyset$. Hence $Z$ is excluded from $G_epsilon$ and Stage 2.
  Because $Delta > 2 epsilon$, no core component intersects both $A$ and
  $B$, so the minimum inter-component distance between any core in $A$
  and any core in $B$ is at least $Delta$. Thus any reducible linkage
  satisfying $D (U, V) >= D_"SL" (U, V)$ (including Single, Average,
  Complete, and Ward) maintains $A$ and $B$ distinct up to merge height
  at least $Delta$, yielding ratio $Delta / delta -> infinity$ as
  $delta -> 0$.
]

#subsec("B", "Reducibility, Ultrametricity, and Permutation Invariance")

#lemma(
  title: "Lance--Williams Reducibility and Permutation Invariance",
  origin: "Original Proof",
)[
  Suppose the Stage-2 Lance--Williams parameters satisfy
  $
    alpha_A + gamma >= 0, quad alpha_B + gamma >= 0, \
    alpha_A + alpha_B >= 1, quad alpha_A + alpha_B + beta >= 1.
  $
  Then: (i) $D$ satisfies Bruynooghe's *reducibility property*
  $D (A, B) <= min(D (A, C), D (B, C)) => D (A union B, C) >=
  min(D (A, C), D (B, C))$, ensuring non-decreasing merge heights
  (ultrametricity on $cal(M)$) and exact $O(m^2)$ nearest-neighbor-chain
  execution @mullner2011modern; and (ii) assuming generic distinct
  border-to-core distances, HCA-DBSCAN is strictly invariant to all
  input permutations of $X$, resolving classical DBSCAN's border-order
  non-determinism @schubert2017dbscanrevisited.
]

#proof[
  (i) Let $u = D (A, C)$, $v = D (B, C)$, and $w = D (A, B)$ satisfy
  $0 <= w <= min(u, v)$. Assume without loss of generality that
  $u <= v$, so $min(u, v) = u$ and $|u - v| = v - u >= 0$. Rewriting
  the Lance--Williams recurrence gives
  $
    D (A union B, C)
      &= (alpha_A - gamma) u + (alpha_B + gamma) v + beta w \
      &>= (alpha_A - gamma) u + (alpha_B + gamma) u + beta w \
      &= (alpha_A + alpha_B) u + beta w,
  $
  where the inequality uses $alpha_B + gamma >= 0$ and $v >= u$. If
  $beta >= 0$, then $beta w >= 0$ yields
  $D (A union B, C) >= (alpha_A + alpha_B) u >= u$. If $beta < 0$, then
  $w <= u$ gives $beta w >= beta u$, whence
  $D (A union B, C) >= (alpha_A + alpha_B + beta) u >= u$. The case
  $v <= u$ is symmetric via $alpha_A + gamma >= 0$. By induction along
  merges, merge heights are non-decreasing (defining an exact
  ultrametric on $cal(M)$) and reciprocal nearest neighbors are
  preserved, validating the $O(m^2)$ NN-chain algorithm
  @mullner2011modern. Note that Single ($alpha = 1/2, beta = 0,
  gamma = -1/2$), Complete ($alpha = 1/2, beta = 0, gamma = 1/2$),
  Average ($alpha_A = abs(A)/(abs(A)+abs(B)), beta = gamma = 0$), and
  Ward ($alpha_A = (abs(A)+abs(C))/(abs(A)+abs(B)+abs(C)), beta =
  -abs(C)/(abs(A)+abs(B)+abs(C)), gamma = 0$) all satisfy these four
  inequalities.

  (ii) In Stage 1, the sets $X_"core"$, $X_"border"$, $X_"noise"$, and
  the connected components $cal(M)$ of $G_epsilon$ depend solely on
  pairwise distances and are symmetric functions of $X$. In Stage 2,
  reducibility guarantees that merging any reciprocal nearest-neighbor
  pair produces the same unique dendrogram (up to tied distances). In
  Stage 3, $x^* (y) = arg min_(x in X_"core") d (y, x)$ depends only on
  the unordered set $X_"core"$ and not on sample ordering. Hence the
  complete clustering output is permutation-invariant.
]

#subsec("C", "Density Level-Set Consistency and Multi-Density Separation")

#theorem(
  title: "Empirical Level-Set Recovery and Multi-Density Resolution",
  origin: "Classical — Cited",
)[
  Let $X_1, dots, X_n$ be i.i.d. draws from a continuous density $f$ on
  $bb(R)^d$, and let $v_d$ be the volume of the unit Euclidean ball.
  Under the bandwidth scaling $epsilon_n -> 0$, $mu_n / log n ->
  infinity$, and $mu_n / (n v_d epsilon_n^d) -> lambda > 0$, Stage 1
  recovers with high probability the connected components of the
  population superlevel set $L_f (lambda) = {x : f (x) >= lambda}$ in
  the sense of Hartigan @hartigan1981consistency and Chaudhuri et al.
  @chaudhuri2014consistent @campello2015hierarchical. Moreover, if $f$
  has two dense modes of peak height $lambda_1$ separated by a valley
  of height $lambda_v < lambda_1$ alongside a sparser cluster of height
  $lambda_2 in (lambda_v, lambda_1)$, setting $lambda in (lambda_v,
  lambda_1)$ isolates both dense modes as separate cores while Stage 2
  re-links fragmented sub-cores of the sparser cluster via hierarchical
  linkage.
]

#proof[
  For any $x in bb(R)^d$, the indicator variables
  $Z_i (x) = bb(I)(norm(X_i - x)_2 <= epsilon_n)$ are i.i.d. Bernoulli
  with success probability $p_n (x) = integral_(B(x, epsilon_n)) f (u)
  dif u = v_d epsilon_n^d (f (x) + o(1))$ uniformly on compact sets
  where $f$ is continuous. Because VC-dimension of Euclidean balls in
  $bb(R)^d$ is $d + 2$, standard uniform Bernstein / VC concentration
  @chaudhuri2014consistent implies that with probability $1 - o(1)$,
  $
    sup_(x in bb(R)^d) |abs(N_(epsilon_n) (x)) / (n v_d epsilon_n^d)
      - f (x)| = O(sqrt((d log n) / (n epsilon_n^d))) + o(1) -> 0.
  $
  Hence $X_"core" = {x : abs(N_(epsilon_n)(x)) / (n v_d epsilon_n^d) >=
  lambda}$ sandwiches the true superlevel set between $L_f (lambda +
  eta_n)$ and $L_f (lambda - eta_n)$ with $eta_n -> 0$. Whenever two
  components of $L_f (lambda)$ are separated by a safety margin
  exceeding $2 epsilon_n$, no edge of $G_(epsilon_n)$ crosses the
  valley, while internal $epsilon_n$-connectivity holds within each
  component @hartigan1981consistency @campello2015hierarchical. In the
  multi-density regime, choosing $epsilon_n$ small enough that
  $lambda > lambda_v$ prevents premature fusion of the two dense modes
  across the valley (where flat DBSCAN tuned to the sparse cluster would
  collapse both dense modes into one), while Stage 2 hierarchically
  merges adjacent sparse-cluster cores before spanning the larger
  inter-cluster gap to the dense modes.
]

#subsec("D", "Exact Huygens--Steiner Variance Compression for Core Ward Linkage")

#theorem(
  title: "Exact Algebraic Equivalence of Core-Compressed Ward Linkage",
  origin: "Original Derivation",
)[
  Let $X subset.eq bb(R)^d$ and let $cal(M) = {M_1, dots, M_m}$ be the
  extracted density cores with weights $n_a = abs(M_a)$ and centroids
  $bar(x)_a = n_a^(-1) sum_(x in M_a) x$. For any macro-partition
  $cal(C) = {C_1, dots, C_K}$ coarsening $cal(M)$,
  $
    "SSE"(X_"core", cal(C)) = sum_(a = 1)^m "SSE"(M_a)
      + sum_(r = 1)^K sum_(M_a subset.eq C_r)
        n_a norm(bar(x)_a - bar(x)_(C_r))_2^2.
  $
  Consequently, Ward linkage computed solely from the $m$ weighted core
  tuples ${(n_a, bar(x)_a)}_(a = 1)^m$ yields the *exact* same objective
  increments $Delta "SSE"(A, B)$ as point-level Ward linkage on
  $X_"core"$ constrained to preserve $cal(M)$, reducing Stage-2 time
  from $O(abs(X_"core")^2 d)$ to $O(m^2 d)$.
]

#proof[
  For any macro-cluster $C_r = union.big_(M_a subset.eq C_r) M_a$ with
  centroid $bar(x)_(C_r)$, expand the squared Euclidean error inside
  each constituent core $M_a$:
  $
    sum_(x in M_a) norm(x - bar(x)_(C_r))_2^2
      &= sum_(x in M_a) norm((x - bar(x)_a) + (bar(x)_a - bar(x)_(C_r)))_2^2 \
      &= "SSE"(M_a) + n_a norm(bar(x)_a - bar(x)_(C_r))_2^2 \
      &quad + 2 chevron.l sum_(x in M_a) (x - bar(x)_a),
        bar(x)_a - bar(x)_(C_r) chevron.r.
  $
  By definition of $bar(x)_a$, $sum_(x in M_a) (x - bar(x)_a) = 0$, so
  the cross-term vanishes identically. Summing over all $M_a subset.eq
  C_r$ and all $r = 1, dots, K$ establishes the decomposition. Because
  $sum_(a = 1)^m "SSE"(M_a)$ is invariant across all coarsenings of
  $cal(M)$, merging two macro-clusters $A, B$ incurs exact Ward cost
  $
    Delta "SSE"(A, B) = (abs(A) abs(B)) / (abs(A) + abs(B))
      norm(bar(x)_A - bar(x)_B)_2^2,
  $
  computable in $O(d)$ from ${(n_a, bar(x)_a)}$ without revisiting raw
  points. Using NN-chain @mullner2011modern on $m$ cores requires
  $O(m^2 d)$ time—a speedup factor of $(abs(X_"core") / m)^2$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

HCA-DBSCAN's runtime splits across Stage 1 (neighborhood queries and
BFS) and Stage 2 (agglomerative linkage on $m$ cores). In general metric
spaces without indexing, Stage 1 costs $O(n^2)$ distance evaluations;
in low-dimensional Euclidean spaces equipped with $k$-d trees, ball
trees, or R-trees, range queries drop to $O(n log n)$ expected time
@ester1996dbscan @schubert2017dbscanrevisited. Stage 2 costs $O(m^2)$
using NN-chain for reducible linkages (Lemma 1) or $O(m^2 log m)$ with
heap-based linkage, plus $O(abs(X_"border") m_"core")$ (or
$O(abs(X_"border") log n)$ via nearest-core index) for Stage 3. When
$m << n$, Stage 2 takes negligible time relative to Stage 1, yielding
$O(n log n + m^2)$ total complexity and $O(n + m^2)$ memory (or $O(n)$
with on-the-fly NN-chain distances).

Geometrically, HCA-DBSCAN exhibits a dual-scale inductive bias: below
scale $epsilon$, clusters follow arbitrary non-convex density-connected
manifolds; above scale $epsilon$, macro-clusters obey the chosen
Stage-2 linkage geometry (e.g., elongated connectivity under Single
Linkage, compact variance balance under core Ward, or moderate bridging
under Average Linkage). Crucially, stripping $X_"noise"$ before Stage 2
insulates the dendrogram from chaining and outlier inflation.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Input Domain], [Metric space $(X, d)$ (Euclidean for Ward Stage 2)],
    [Time Complexity], [$O(n^2)$ metric / $O(n log n + m^2)$ spatial-indexed],
    [Space Complexity], [$O(n + m^2)$ (or $O(n)$ via on-the-fly NN-chain)],
    [Cluster Geometry], [Non-convex density cores + linkage macro-structure],
    [Key Hyperparameters], [Core radius $epsilon$, density $mu$, Stage-2 linkage],
    [Noise Handling], [Explicit noise rejection ($X_"noise"$) prior to linkage],
    [Ultrametricity], [Exact on $cal(M)$ for reducible Stage-2 linkages],
    [Invariances], [Strict permutation invariance; isometric invariance],
  ),
  caption: [Operational and theoretical profile of HCA-DBSCAN.],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

HCA-DBSCAN is particularly effective whenever datasets combine
background clutter, irregular local manifolds, and hierarchical
multi-density macro-structure:

+ *LiDAR Point-Cloud and Autonomous Driving Perception:* 3D automotive
  LiDAR scans contain dense surfaces (vehicles, pedestrians, curbs)
  surrounded by atmospheric scatter and ground clutter. Stage 1 strips
  sparse scatter into $X_"noise"$ and compresses dense object patches
  into cores $cal(M)$, while Stage 2 merges multi-part objects (e.g.,
  partially occluded vehicles) hierarchically without chaining across
  road noise.
+ *Geospatial Mobility and Multi-Scale Hotspot Detection:* Urban GPS
  check-ins form dense downtown cores alongside sparser suburban hubs.
  Single-scale DBSCAN either fuses adjacent downtown districts or drops
  suburbs as noise; HCA-DBSCAN extracts fine cores with tight $epsilon$
  and reconstructs neighborhood-to-district hierarchies in Stage 2.
+ *Customer Segmentation with Outlier Transactions:* Retail behavioral
  tables mix tightly clustered loyal buyer segments with erratic
  one-off or fraudulent accounts. Purging erratic accounts in Stage 1
  stabilizes Stage-2 Average or Ward linkage, preventing extreme spenders
  from distorting customer archetypes.
+ *Single-Molecule Localization Microscopy (SMLM):* Super-resolution
  fluorescence microscopy localizes protein nanoclusters amidst random
  background blinking artifacts. HCA-DBSCAN rejects isolated background
  localizations via $(epsilon, mu)$ filtering and builds biologically
  interpretable multi-scale protein-complex hierarchies across core
  aggregates.
+ *Cybersecurity Flow and Log-Event Triage:* Network intrusion flows
  cluster into dense attack campaigns embedded in noisy background
  traffic. Core compression shrinks millions of flows to $m$ compact
  campaign signatures, enabling rapid hierarchical campaign attribution
  via Stage-2 linkage.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four operational failure regimes bound HCA-DBSCAN. First, *overly large
$epsilon$ collapse*: if $epsilon$ exceeds the narrowest inter-cluster
gap, Stage 1 merges distinct modes into a single core component before
Stage 2 can act, rendering Stage 2 powerless to split them. Second,
*overly small $epsilon$ or large $mu$ starvation*: if sparse true
clusters fail the core threshold $abs(N_epsilon (x)) >= mu$, they are
discarded permanently into $X_"noise"$. Third, *high-dimensional curse
of dimensionality* @beyer1999nearest: volume $v_d epsilon^d$ scales
exponentially with $d$, making fixed-$epsilon$ ball density brittle and
degrading spatial-index acceleration back to $O(n^2)$. Fourth, *Stage-2
linkage mismatch*: applying Ward linkage to deeply intertwined spirals
whose cores were fragmented at small $epsilon$ can mis-merge opposite
arms by centroid proximity; Single or Average Stage-2 linkage is needed
for non-convex macro-shapes.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Non-convex + dense noise], [Ideal: noise purged before core linkage], [Keep $epsilon$ conservative; use SL/Avg Stage 2],
    [Multi-scale customer hubs], [Strong: core compression + Ward/Avg hierarchy], [Check $X_"noise"$ fraction for rare tiers],
    [Extreme density contrast], [Moderate: very sparse modes risk noise drop], [Use HDBSCAN @campello2013hdbscan or OPTICS @ankerst1999optics],
    [High-D ($d > 30$) vectors], [Weak: $epsilon$-ball density concentrates], [Apply PCA/UMAP or subspace methods],
    [Zero-noise spherical data], [Redundant Stage 1; slight border distortion], [Use Ward, BUBBLE, or $k$-means],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for HCA-DBSCAN.],
)

#sec("V", "Structural Relations and Empirical Behavior")

HCA-DBSCAN belongs to the family of *two-phase density-hierarchical
hybrids* alongside CURE @guha1998cure, Chameleon @karypis1999chameleon,
OPTICS @ankerst1999optics, and HDBSCAN @campello2013hdbscan
@campello2015hierarchical. Whereas HDBSCAN eliminates the global
parameter $epsilon$ by replacing raw distances with mutual reachability
distances $d_"mreach" (x, y) = max(d_"core" (x), d_"core" (y),
d (x, y))$ and extracting excess-of-mass stable clusters, HCA-DBSCAN
fixes an explicit density resolution $epsilon$ to achieve hard noise
pruning and core compression ($n -> m$), then allows arbitrary
Lance--Williams Stage-2 criteria (including Ward and Average linkage)
that HDBSCAN's single-linkage mutual-reachability tree cannot express.
Similarly, DBSCAN++ @jang2019dbscanpp accelerates Stage-1 core discovery
via subsampled cores, complementing HCA-DBSCAN's Stage-2 agglomeration.

#proposition(
  title: "Relation Between HCA-DBSCAN and HDBSCAN Mutual Reachability",
  origin: "Original Proof",
)[
  Let $d_("core", mu) (x)$ denote the distance from $x$ to its $mu$-th
  nearest neighbor in $X$ @campello2013hdbscan. Then $X_"core" =
  {x in X : d_("core", mu) (x) <= epsilon}$, and for every pair of core
  points $x, y in X_"core"$ belonging to distinct density cores
  $M_a != M_b$, the HDBSCAN mutual reachability distance coincides
  *identically* with the raw metric distance:
  $
    d_("mreach", mu) (x, y) = max(d_("core", mu) (x), d_("core", mu) (y),
      d (x, y)) = d (x, y).
  $
]

#proof[
  By definition, $abs(N_epsilon (x)) >= mu$ iff at least $mu$ points
  (including $x$ itself) lie within distance $epsilon$ of $x$, which is
  equivalent to $d_("core", mu) (x) <= epsilon$. If $x in M_a$ and
  $y in M_b$ with $a != b$, then $(x, y) in.not E_epsilon$, so
  $d (x, y) > epsilon >= max(d_("core", mu) (x), d_("core", mu) (y))$.
  Hence $max(d_("core", mu) (x), d_("core", mu) (y), d (x, y)) = d (x, y)$.
]

Proposition 1 reveals an exact mathematical bridge: on the extracted
core set $X_"core"$, HCA-DBSCAN with Stage-2 Single Linkage agrees
edge-for-edge with HDBSCAN's mutual-reachability MST above scale
$epsilon$, while cleanly excising all points with $d_("core", mu) (x) >
epsilon$ prior to hierarchy construction.

The ClusterOverview 4-dataset benchmark strikingly validates this
two-stage design @rousseeuw1987silhouettes @hubert1985comparing.
HCA-DBSCAN achieves the strongest overall profile across non-convex,
heterogeneous, and high-dimensional tasks among Algorithms 68--71: on
Multi-Shape it reaches an exceptional Silhouette of *0.804*
($"ARI" = "N/A"$)—surpassing Genie (0.718), BUBBLE (0.452), and MONA
(0.304)—directly reflecting Theorem 1's noise-bridge immunity on
irregular geometries. On Customers, HCA-DBSCAN attains $"Sil" = 0.514$
(dramatically exceeding Genie's 0.266 and BUBBLE's 0.353), because
Stage 1 isolates heavy-tail customer outliers before Stage 2 merges
unequal-size behavioral cores without artificial cardinality pressure.
On High-D it leads the group at $"Sil" = 0.324$ ($"ARI" = "N/A"$), while
on Iris it records $"Sil" = 0.336$ and $"ARI" = 0.545$ (slightly below
Genie's 0.606 ARI due to overlapping versicolor/virginica density
borders). Finally, core compression ($m << n$) drives a lean mean
runtime of *1.008 ms*, nearly $2.8 times$ faster than Genie (2.781 ms)
and faster than BUBBLE (1.215 ms).
