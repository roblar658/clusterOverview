#import "../template/ieee.typ": *

#alg_entry(
  id: 65,
  name: "CURE (Clustering Using REpresentatives)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1998",
  paradigm: "Shrunk Multi-Representative Agglomerative Clustering",
  complexity: $cal(O)(n_s^2 log n_s + n k c d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Classical agglomerative clustering sits between two brittle geometric extremes: *centroid-based* algorithms (e.g., $k$-means, Ward @ward1963hierarchical, and coarse-threshold BIRCH @zhang1996birch) represent each cluster by a single mean point and therefore fail on elongated, non-spherical, or widely varying-size clusters; conversely, *all-points* single-linkage clustering @sibson1973slink captures arbitrary shapes but collapses under a thin bridge of noise points (the *chaining effect*). CURE @guha1998cure resolves this dichotomy by representing each cluster with a constant number $c >= 1$ of well-scattered boundary points shrunk uniformly toward the cluster centroid by a factor $alpha in [0, 1]$, paired with Chernoff-bounded random sampling and partitioned pre-clustering for large datasets $X = { bold(x)_1, dots, bold(x)_n } subset bb(R)^d$.

#definition(title: "Shrunk Representative Set and CURE Dissimilarity")[
  Let $C subset bb(R)^d$ be a nonempty finite cluster with centroid $bold(mu)(C) = (1 / abs(C)) sum_(bold(x) in C) bold(x)$, and let $c in bb(N)_+$ and $alpha in [0, 1]$. The *well-scattered seed set* $S_c (C) = { bold(p)_1, dots, bold(p)_tilde(c) } subset.eq C$ (with $tilde(c) = min(c, abs(C))$) is generated greedily via farthest-first traversal:
  $ bold(p)_1 &= arg max_(bold(x) in C) norm(bold(x) - bold(mu)(C))_2, \
    bold(p)_r &= arg max_(bold(x) in C) min_(1 <= j < r) norm(bold(x) - bold(p)_j)_2, quad r = 2, dots, tilde(c). $
  The *$alpha$-shrunk representative set* of $C$ is the affine contraction toward $bold(mu)(C)$:
  $ R_alpha (C) = { bold(q)_r = (1 - alpha) bold(p)_r + alpha bold(mu)(C) : bold(p)_r in S_c (C) }. $
  The *CURE inter-cluster dissimilarity* between clusters $U, V subset bb(R)^d$ is the minimum Euclidean distance between their shrunk representatives:
  $ d_"CURE" (U, V; c, alpha) = min_(bold(q) in R_alpha (U), bold(q)' in R_alpha (V)) norm(bold(q) - bold(q)')_2. $
]

Geometrically, $S_c (C)$ captures the extremal geometry and orientation of $C$ (e.g., endpoints and ridges of elongated ellipsoids), while shrinking toward $bold(mu)(C)$ by $alpha$ pulls isolated outlier points radially inward much further than core boundary points, dampening chaining without collapsing cluster geometry to a single point.

#algbox(
  title: "CURE Scalable Multi-Representative Agglomeration",
  inputs: [Dataset $X subset bb(R)^d$ ($abs(X) = n$), target clusters $k$, representative count $c$, shrinkage $alpha in [0, 1]$, sample size $n_s$, partition count $p$, reduction factor $q > 1$.],
  outputs: [Final representative sets ${ R_alpha (C_1), dots, R_alpha (C_k) }$ and full-dataset assignment $pi: X -> { 1, dots, k }$.],
)[
  + *Random Sampling:* Draw a uniform random sample $S subset.eq X$ of size $n_s$ sized via the Chernoff bound (Theorem 3).
  + *Partitioned Pre-Clustering:* Split $S$ into $p$ equal partitions of size $n_s / p$. Run the core agglomeration loop (Step 4) independently on each partition until $n_s / (p q)$ clusters remain per partition, then discard sparse outlier clusters.
  + *Index Initialization:* Pool all surviving pre-clusters into a single list. Insert every representative $bold(q) in R_alpha (u)$ into a $k$-d tree $cal(T)$, and insert each cluster $u$ into a min-heap $Q$ keyed by $d_"CURE" (u, u."closest")$.
  + *Core Agglomeration Loop:* While $abs(Q) > k$:
    - Extract $u = arg min_(v in Q) d_"CURE" (v, v."closest")$ and let $v = u."closest"$.
    - Form $w = u union v$, compute $bold(mu)(w) = (abs(u) bold(mu)(u) + abs(v) bold(mu)(v)) / (abs(u) + abs(v))$, select $S_c (w)$ via farthest-first traversal on $S_c (u) union S_c (v)$ (or $w$), and construct $R_alpha (w)$.
    - Delete representatives of $u, v$ from $cal(T)$, insert $R_alpha (w)$ into $cal(T)$, and update `closest` pointers and heap keys in $Q$ via $k$-d tree range/NN queries.
    - At an intermediate phase ($abs(Q) approx k_"prune"$), eliminate slow-growing outlier clusters of tiny cardinality.
  + *Full-Dataset Labeling:* Stream all remaining points $bold(x) in X without S$ from disk and assign $pi(bold(x)) = arg min_(j in { 1, dots, k }) min_(bold(q) in R_alpha (C_j)) norm(bold(x) - bold(q))_2$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Continuous Homotopy Between Single-Linkage and Centroid-Linkage")

We first prove that CURE's dissimilarity defines a Lipschitz-continuous, similarity-equivariant homotopy interpolating strictly between single-linkage and centroid-linkage.

#theorem(title: "Boundary Reductions, Lipschitz Continuity, and Equivariance", origin: "Original Proof")[
  For any nonempty finite clusters $U, V subset bb(R)^d$ and integer $c >= 1$, let $R_"max" (C) = max_(bold(x) in C) norm(bold(x) - bold(mu)(C))_2$. Then:
  + *Centroid Limit ($alpha = 1$):* For every $c >= 1$, $d_"CURE" (U, V; c, 1) = norm(bold(mu)(U) - bold(mu)(V))_2$.
  + *Single-Linkage Limit ($alpha = 0$):* Whenever $c >= max(abs(U), abs(V))$, $d_"CURE" (U, V; c, 0) = min_(bold(x) in U, bold(y) in V) norm(bold(x) - bold(y))_2 = d_"SL" (U, V)$.
  + *Lipschitz Continuity in $alpha$:* For any $alpha_1, alpha_2 in [0, 1]$,
    $ &abs(d_"CURE" (U, V; c, alpha_1) - d_"CURE" (U, V; c, alpha_2)) \
      &<= abs(alpha_1 - alpha_2) (R_"max" (U) + R_"max" (V)). $
  + *Similarity Equivariance:* Under any affine similarity $Phi(bold(x)) = s Q bold(x) + bold(b)$ ($s > 0$, $Q^top Q = I_d$, $bold(b) in bb(R)^d$), $d_"CURE" (Phi(U), Phi(V)) = s dot d_"CURE" (U, V)$, making the unpruned merge dendrogram invariant to translation, rotation, reflection, and uniform scaling.
]

#proof[
  For Part (1), setting $alpha = 1$ collapses every representative $bold(q)_r = (1 - 1) bold(p)_r + 1 dot bold(mu)(C) = bold(mu)(C)$, so $R_1 (U) = { bold(mu)(U) }$ and $R_1 (V) = { bold(mu)(V) }$, giving $d_"CURE" = norm(bold(mu)(U) - bold(mu)(V))_2$.

  For Part (2), when $c >= max(abs(U), abs(V))$, farthest-first traversal exhausts all elements, so $S_c (U) = U$ and $S_c (V) = V$. Setting $alpha = 0$ yields $R_0 (U) = U$ and $R_0 (V) = V$, reducing $d_"CURE"$ identically to single-linkage $d_"SL" (U, V)$.

  For Part (3), note that $S_c (U)$ and $S_c (V)$ are independent of $alpha$. For any fixed pair $(bold(p), bold(p)') in S_c (U) times S_c (V)$, denote $bold(q)(alpha) = (1 - alpha) bold(p) + alpha bold(mu)(U)$ and $bold(q)'(alpha) = (1 - alpha) bold(p)' + alpha bold(mu)(V)$. By the reverse triangle inequality,
  $ &abs(norm(bold(q)(alpha_1) - bold(q)'(alpha_1))_2 - norm(bold(q)(alpha_2) - bold(q)'(alpha_2))_2) \
    &<= norm((bold(q)(alpha_1) - bold(q)(alpha_2)) - (bold(q)'(alpha_1) - bold(q)'(alpha_2)))_2 \
    &= abs(alpha_1 - alpha_2) norm((bold(mu)(U) - bold(p)) - (bold(mu)(V) - bold(p)'))_2 \
    &<= abs(alpha_1 - alpha_2) (R_"max" (U) + R_"max" (V)). $
  Using the elementary bound $abs(min_i a_i - min_i b_i) <= max_i abs(a_i - b_i)$ over all $(bold(p), bold(p)') in S_c (U) times S_c (V)$ proves Part (3).

  For Part (4), $Phi$ preserves centroids ($bold(mu)(Phi(C)) = Phi(bold(mu)(C))$) and scales all pairwise distances uniformly by $s > 0$, preserving every argmax in farthest-first traversal (up to deterministic tie-breaking). Moreover,
  $ (1 - alpha) Phi(bold(p)) + alpha Phi(bold(mu)(C)) = Phi((1 - alpha) bold(p) + alpha bold(mu)(C)), $
  so $R_alpha (Phi(C)) = Phi(R_alpha (C))$ and all inter-representative distances scale identically by $s$, preserving the complete merge order.
]

#subsec("B", "Quantitative Anti-Chaining Margin Under Outlier Shrinkage")

A principal heuristic claim of CURE is that shrinking representatives toward the centroid suppresses single-linkage chaining caused by bridge outliers. We make this guarantee mathematically exact.

#theorem(title: "Outlier Shrinkage Margin and Chaining Suppression", origin: "Original Proof")[
  Let $U, V subset bb(R)^d$ be compact clusters contained in Euclidean balls $B_2 (bold(mu)_U, r_U)$ and $B_2 (bold(mu)_V, r_V)$ with centroid separation $Delta_mu = norm(bold(mu)_U - bold(mu)_V)_2 > r_U + r_V$ and cardinality $n_U = abs(U) >= 1$. Suppose an adversarial bridge outlier $bold(z) in bb(R)^d$ at distance $L = norm(bold(z) - bold(mu)_U)_2 >= r_U$ is absorbed into $U' = U union { bold(z) }$. Then for every $c >= 1$ and $alpha in [0, 1]$,
  $ d_"CURE" (U', V; c, alpha) >= Delta_mu - (1 - alpha) r_V - (1 - alpha n_U / (n_U + 1)) L. $
  In particular, even when $bold(z)$ touches the boundary of $V$ ($L = Delta_mu - r_V$, driving single-linkage separation to $d_"SL" (U', V) = 0$), CURE guarantees a strictly positive anti-chaining margin
  $ d_"CURE" (U', V; c, alpha) >= alpha (n_U / (n_U + 1) L + r_V) > 0 quad forall alpha in (0, 1]. $
]

#proof[
  The perturbed centroid of $U' = U union { bold(z) }$ is
  $ bold(mu)_(U') = (n_U bold(mu)_U + bold(z)) / (n_U + 1) = bold(mu)_U + 1 / (n_U + 1) (bold(z) - bold(mu)_U). $
  We bound the radial displacement $norm(bold(q) - bold(mu)_U)_2$ of every shrunk representative $bold(q) in R_alpha (U')$:
  - *Case 1 ($bold(p) = bold(z)$):* The shrunk outlier satisfies
    $ bold(q)_z - bold(mu)_U &= (1 - alpha)(bold(z) - bold(mu)_U) + alpha (bold(mu)_(U') - bold(mu)_U) \
      &= (1 - alpha + alpha / (n_U + 1)) (bold(z) - bold(mu)_U) \
      &= (1 - alpha n_U / (n_U + 1)) (bold(z) - bold(mu)_U), $
    whose Euclidean norm equals $(1 - alpha n_U / (n_U + 1)) L$.
  - *Case 2 ($bold(p) in U$):* Since $norm(bold(p) - bold(mu)_U)_2 <= r_U <= L$, the triangle inequality gives
    $ norm(bold(q)_p - bold(mu)_U)_2 &<= (1 - alpha) norm(bold(p) - bold(mu)_U)_2 + alpha norm(bold(mu)_(U') - bold(mu)_U)_2 \
      &<= (1 - alpha) r_U + alpha / (n_U + 1) L \
      &<= (1 - alpha + alpha / (n_U + 1)) L = (1 - alpha n_U / (n_U + 1)) L. $
  Hence every $bold(q) in R_alpha (U')$ satisfies $norm(bold(q) - bold(mu)_U)_2 <= (1 - alpha n_U / (n_U + 1)) L$. Simultaneously, every $bold(q)' in R_alpha (V)$ satisfies $norm(bold(q)' - bold(mu)_V)_2 = (1 - alpha) norm(bold(p)' - bold(mu)_V)_2 <= (1 - alpha) r_V$. Applying the reverse triangle inequality to $bold(q) - bold(q)' = (bold(mu)_U - bold(mu)_V) + (bold(q) - bold(mu)_U) - (bold(q)' - bold(mu)_V)$ yields
  $ norm(bold(q) - bold(q)')_2 >= Delta_mu - (1 - alpha n_U / (n_U + 1)) L - (1 - alpha) r_V. $
  Taking the minimum over $R_alpha (U') times R_alpha (V)$ gives the general lower bound. Substituting $Delta_mu = L + r_V$ simplifies the right-hand side directly to $alpha ((n_U / (n_U + 1)) L + r_V) > 0$.
]

Theorem 2 reveals why CURE's shrinkage factor operates asymmetrically: the radial pull-back on a distant outlier is $alpha (n_U / (n_U + 1)) L$, proportional to $L$ itself, while core points retract only by $alpha r_U$.

#subsec("C", "Chernoff Sample-Size Bound for Cluster Preservation")

To avoid $Omega(n^2)$ runtime on massive datasets, CURE clusters a uniform random sample $S subset.eq X$ of size $n_s$. The following classical bound certifies the minimum $n_s$ guaranteeing that every cluster is faithfully represented.

#theorem(title: "Chernoff Sampling Bound for Cluster Representation", origin: "Classical — Cited")[
  @guha1998cure Let $X$ of size $n$ contain $k$ clusters $C_1, dots, C_k$ with minimum cluster size $n_"min" = min_(1 <= j <= k) abs(C_j)$, and fix a target retention fraction $f in (0, 1)$ and failure probability $delta in (0, 1)$. If a uniform random sample $S$ of size $n_s$ (drawn with or without replacement) satisfies
  $ n_s &>= f n + n / n_"min" ln(k / delta) \
    &quad + n / n_"min" sqrt((ln(k / delta))^2 + 2 f n_"min" ln(k / delta)), $
  then with probability at least $1 - delta$, *every* cluster simultaneously satisfies $abs(S inter C_j) >= f abs(C_j)$.
]

#proof[
  Fix a cluster $C_j$ of size $n_j >= n_"min"$ and set $p_j = n_j / n$. Under sampling with replacement, $Y_j = abs(S inter C_j) tilde "Bin"(n_s, p_j)$ with mean $mu_j = n_s p_j$. (Under sampling without replacement, $Y_j$ is hypergeometric, which is convex-order dominated by $"Bin"(n_s, p_j)$, so the identical Chernoff upper tail holds.) Assuming $n_s > f n$, let $gamma = 1 - f n_j / mu_j = 1 - f n / n_s in (0, 1)$. By the multiplicative lower-tail Chernoff inequality,
  $ "Pr"(Y_j < f n_j) = "Pr"(Y_j < (1 - gamma) mu_j) <= exp(- (gamma^2 mu_j) / 2). $
  Requiring $exp(-(gamma^2 mu_j) / 2) <= delta / k$ is algebraically equivalent to
  $ (n_s - f n)^2 / (2 n_s) dot n_j / n >= ln(k / delta). $
  Let $a = f n$ and $b_j = (n / n_j) ln(k / delta)$. Expanding gives the quadratic inequality $n_s^2 - 2 (a + b_j) n_s + a^2 >= 0$, whose upper root (the unique root exceeding $a = f n$) is
  $ n_s >= a + b_j + sqrt(b_j^2 + 2 a b_j). $
  Since $n_j >= n_"min"$ implies $b_j <= b_"min" = (n / n_"min") ln(k / delta)$, substituting $b_"min"$ ensures $"Pr"(Y_j < f n_j) <= delta / k$ for every $j$. Applying a union bound across $j = 1, dots, k$ gives $"Pr"(exists j: Y_j < f n_j) <= sum_(j=1)^k delta / k = delta$.
]

Notice that the sample size scales inversely with the *smallest* cluster proportion $n_"min" / n$, making $n_s$ independent of total dataset size $n$ once cluster proportions $n_j / n$ are fixed (setting $f$ proportional to $c / n_"min"$ yields $n_s = cal(O)((n / n_"min") (c + ln(k / delta)))$).

#subsec("D", "Concave Geometry Collapse Under Excessive Shrinkage")

While CURE handles elongated convex shapes (ellipsoids, parallel bars), affine shrinkage toward the *Euclidean centroid* implicitly assumes each cluster is star-shaped around its mean. For non-convex or concentric manifolds whose centroid falls inside a neighboring cluster, shrinkage causes provable geometric inversion.

#proposition(title: "Critical Shrinkage Collapse on Concentric Manifolds", origin: "Original Counterexample")[
  Consider two clusters in $bb(R)^2$: an outer ring $U = { bold(x) in bb(R)^2 : norm(bold(x))_2 = R }$ and a concentric inner disk $V = { bold(x) in bb(R)^2 : norm(bold(x))_2 <= r }$ with $0 < r < R$, both having centroid $bold(mu)_U = bold(mu)_V = bold(0)$.
  + *Inter-Cluster Contraction:* $d_"CURE" (U, V; c, alpha) = (1 - alpha)(R - r)$, collapsing linearly to $0$ as $alpha -> 1$.
  + *Phase 5 Misclassification Threshold:* For every shrinkage factor exceeding the critical threshold
    $ alpha > alpha^* = (R - r) / (R + r), $
    every boundary point $bold(y) in partial V$ aligned radially with a representative of $U$ is strictly closer to $R_alpha (U)$ than to $R_alpha (V)$, causing Phase 5 nearest-representative labeling to misclassify the boundary of $V$ into $U$.
]

#proof[
  By rotational symmetry, $bold(mu)_U = bold(mu)_V = bold(0)$. Every scattered point $bold(p) in S_c (U)$ has $norm(bold(p))_2 = R$, so its shrunk representative $bold(q) = (1 - alpha) bold(p)$ lies on the circle of radius $(1 - alpha) R$. Every $bold(p)' in S_c (V)$ selected by farthest-first traversal (starting from $partial V$) lies at radius $r$, so $norm(bold(q)')_2 = (1 - alpha) r$. For aligned points, the minimum distance between circles of radii $(1 - alpha) R$ and $(1 - alpha) r$ is $(1 - alpha)(R - r)$, proving Part (1).

  For Part (2), consider a boundary point $bold(y) = r hat(bold(u)) in partial V$ aligned with unit direction $hat(bold(u))$ of representatives $bold(q) = (1 - alpha) R hat(bold(u)) in R_alpha (U)$ and $bold(q)' = (1 - alpha) r hat(bold(u)) in R_alpha (V)$ (note that any other representative of $V$ lies at distance $>= norm(bold(y) - bold(q)')_2$). The distances are
  $ norm(bold(y) - bold(q)')_2 &= r - (1 - alpha) r = alpha r, \
    norm(bold(y) - bold(q))_2 &= abs((1 - alpha) R - r). $
  We check when $norm(bold(y) - bold(q))_2 < norm(bold(y) - bold(q)')_2$:
  - If $(1 - alpha) R >= r$, this reduces to $(1 - alpha) R - r < alpha r <=> (1 - alpha) R < (1 + alpha) r <=> alpha > (R - r) / (R + r)$.
  - If $(1 - alpha) R < r$ (i.e., $alpha > 1 - r / R > (R - r) / (R + r)$), then $r - (1 - alpha) R < alpha r <=> (1 - alpha) r < (1 - alpha) R$, which holds identically for all $alpha < 1$ since $r < R$.
  Thus for all $alpha in ((R - r) / (R + r), 1)$, $bold(y) in partial V$ is strictly closer to $R_alpha (U)$ than to $R_alpha (V)$ and is misassigned to $U$.
]

Even at a modest radius ratio $r / R = 0.7$, Proposition 1 yields $alpha^* = 3 / 17 approx 0.176$, proving that CURE's canonical default $alpha = 0.2$ already inverts tightly nested concentric rings.

#sec("III", "Complexity, Geometry and Operational Profile")

*Time Complexity:* On a sample of size $n_s$, building the initial $k$-d tree takes $cal(O)(n_s d log n_s)$. In low dimensions ($d <= 4$), each $k$-d tree nearest-neighbor search averages $cal(O)(log n_s)$, and over $n_s - k$ merges—each updating $cal(O)(c)$ representatives and heap entries—CURE runs in $cal(O)(n_s^2 log n_s)$ worst-case time @guha1998cure (because each merge can invalidate up to $cal(O)(n_s)$ clusters whose `closest` pointer was the merged pair, requiring $cal(O)(n_s)$ heap updates). Partitioning $S$ into $p$ parts of size $n_s / p$ accelerates pre-clustering by a factor of $approx p$. Phase 5 streams the remaining $n - n_s$ points against $k c$ representatives in $cal(O)(n k c d)$ time. Crucially, in high dimensions ($d >> 10$), $k$-d tree backtracking degenerates to linear scan @beyer1999nearest, inflating sample clustering toward $cal(O)(n_s^2 c d)$.

*Space Complexity:* CURE stores only the sample $S$, $k$-d tree, and min-heap in RAM, requiring $cal(O)(n_s c d)$ memory—strictly independent of full dataset size $n$ once $n_s$ is fixed via Theorem 3.

*Geometry and Tuning:* Setting $c approx 10$ and $alpha in [0.2, 0.3]$ captures arbitrary convex ellipsoids, unequal-size clusters, and mildly curved ribbons while suppressing single-linkage chaining (Theorem 2). However, unlike BIRCH's exact sufficient statistics, CURE's farthest-first representative selection is sensitive to the sample quality and cannot be updated purely algebraically from low-order moments.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Sample Phase Time], [$cal(O)(n_s^2 log n_s)$ (low $d$); $cal(O)(n_s^2 c d)$ (high-$d$ $k$-d fallback)],
    [Full Labeling Time], [$cal(O)(n k c d)$ single sequential disk pass],
    [Memory Footprint], [$cal(O)(n_s c d)$ in RAM, independent of $n$],
    [Inductive Bias], [Star-shaped / elongated clusters with unequal sizes],
    [Key Hyperparameters], [Rep. count $c$, shrinkage $alpha$, sample $n_s$, partitions $p, q$, $k$],
    [Invariances], [Translation, orthogonal rotation, and uniform scale invariant],
    [Outlier Handling], [Two-stage sparse-cluster pruning + $alpha$-shrinkage damping],
    [Dendrogram Order], [Not Lance-Williams; can exhibit mild centroid-drift inversions],
  ),
  caption: [Operational, complexity, and geometric profile of CURE (Alg. 65).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Geospatial GIS and Urban Footprint Analysis:* Geographic point clouds (e.g., road corridors, riverine settlements, coastal urban sprawl, $d = 2$) consist of highly elongated, anisotropic, unequal-area clusters surrounded by isolated rural outliers. CURE's $c$ scattered representatives trace linear corridors naturally while $alpha$-shrinkage (Theorem 2) prevents spurious rural bridge mergers.
+ *Astronomical Filament and Galaxy-Group Detection:* 2D/3D sky-survey catalogs contain massive spherical clusters adjacent to small compact dwarf groups and filamentary streams. Unlike centroid methods that split large clusters to absorb adjacent small ones, CURE's multi-representative boundary adapts to extreme size ratios ($100 : 1$).
+ *Customer Behavioral Segmentation with Heavy Outliers:* Retail and telecom datasets often feature imbalanced customer tiers (a large mass-market body next to narrow high-value niches) contaminated by noisy accounts. CURE's Chernoff sampling (Theorem 3) and two-phase sparse-cluster pruning isolate clean segments without letting outliers distort centroids.
+ *Flow Cytometry and Cell-Population Gating:* Low-dimensional ($d = 4$–$12$) fluorescence cytograms exhibit skewed, crescent-like, and elliptical cell populations of widely varying abundances. Multi-representative boundaries conform to skewed biological manifolds far more accurately than spherical $k$-means or BIRCH leaves.
+ *Large-Scale Database Physical Clustering and Indexing:* In spatial databases ($n > 10^7$), CURE's partitioned pre-clustering + single-pass Phase 5 labeling groups spatial objects onto disk pages with minimal random I/O, preserving locality along irregular geographic shapes.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CURE breaks down in four regimes: (i) *Concentric or deeply concave manifolds*, where centroid shrinkage crosses cluster boundaries past $alpha^* = (R - r) / (R + r)$ (Proposition 1); (ii) *High-dimensional feature spaces*, where both $k$-d tree pruning and farthest-first boundary selection degrade under distance concentration @beyer1999nearest; (iii) *Heterogeneous local densities*, where fixed $alpha$ lacks CHAMELEON's relative density normalization; and (iv) *Tiny minority clusters*, where Theorem 3 demands $n_s = Omega(n / n_"min")$, eroding sampling speedup when rare clusters have $n_"min" << n$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [2D/3D geospatial corridors], [Multi-rep tracks elongation; $alpha$ blocks noise], [Keep $alpha in [0.2, 0.3]$, $c >= 10$],
    [Unequal-size convex clusters], [Boundary reps avoid big-cluster splitting], [Ensure $n_s$ covers $n_"min"$ (Thm. 3)],
    [Nested rings / spirals], [Fails above $alpha^* = (R - r)/(R + r)$ (Prop. 1)], [Use CHAMELEON @karypis1999chameleon or HDBSCAN @campello2013hdbscan],
    [High-$d$ vectors ($d > 30$)], [$k$-d tree degrades; distances concentrate], [Use ROCK @guha2000rock or Spectral / SNN @ertoz2003snn],
    [Categorical / set data], [Centroid shrinkage undefined on discrete sets], [Use ROCK @guha2000rock or $k$-Modes @huang1998extensions],
  ),
  caption: [Domain suitability and failure-mode matrix for CURE.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CURE occupies a central bridge in Cat. III hierarchical clustering, unifying classical Lance-Williams linkages @lance1967general with modern scalable sampling while exposing the limitations that motivated graph-based dynamic merging.

#proposition(title: "Linkage Family Embedding and Non-Monotonicity Contrast", origin: "Original Proof")[
  + *Classical Linkage Containment:* Varying $(c, alpha) in bb(N)_+ times [0, 1]$ continuously embeds both single-linkage ($c = n, alpha = 0$ @sibson1973slink) and centroid-linkage ($alpha = 1$ @mullner2011modern) as boundary cases of CURE (Theorem 1), while remaining strictly distinct from complete-linkage @defays1977clink and average-linkage @sokal1958statistical because $d_"CURE"$ minimizes over representatives rather than maximizing or averaging.
  + *Non-Lance-Williams Structure:* For $c >= 2$ and $alpha in (0, 1)$, $d_"CURE" (U union V, W)$ is *not* a function solely of $(d_"CURE" (U, W), d_"CURE" (V, W), d_"CURE" (U, V), abs(U), abs(V), abs(W))$; hence CURE does not satisfy the Lance-Williams recurrence @lance1967general and requires explicit geometric recomputation of $R_alpha (U union V)$.
]

#proof[
  Part (1) follows from Theorem 1(1)–(2); adding an interior point to $U$ never alters $d_"CURE"$ once $S_c (U)$ is fixed, whereas it alters average linkage (and adding a far exterior point alters complete linkage regardless of $alpha$).

  For Part (2), fix $c = 2, alpha = 1 / 2$ on $bb(R)^1$ with $V = { 2 }$ and $W = { 10 }$, and compare two distinct 2-point clusters $U_1 = { -3, -1 }$ and $U_2 = { -6, 0 }$.
  - *Before Merge:* $bold(mu)(U_1) = -2$ yields $R_(1 / 2) (U_1) = { -2.5, -1.5 }$, while $bold(mu)(U_2) = -3$ yields $R_(1 / 2) (U_2) = { -4.5, -1.5 }$. Thus all six Lance-Williams inputs coincide identically: $abs(U_1) = abs(U_2) = 2$, $abs(V) = abs(W) = 1$, $d_"CURE" (U_1, V) = d_"CURE" (U_2, V) = 3.5$, $d_"CURE" (U_1, W) = d_"CURE" (U_2, W) = 11.5$, and $d_"CURE" (V, W) = 8$.
  - *After Merging $U union V$:* For $U_1 union V = { -3, -1, 2 }$, $bold(mu) = -2 / 3$ and $S_2 = { -3, 2 }$, placing the rightmost shrunk representative at $0.5 (2) + 0.5 (-2 / 3) = 2 / 3$, so $d_"CURE" (U_1 union V, W) = 28 / 3$. For $U_2 union V = { -6, 0, 2 }$, $bold(mu) = -4 / 3$ and $S_2 = { -6, 2 }$, placing the rightmost shrunk representative at $0.5 (2) + 0.5 (-4 / 3) = 1 / 3$, so $d_"CURE" (U_2 union V, W) = 29 / 3 != 28 / 3$.
  Hence $d_"CURE"$ cannot satisfy any Lance-Williams recurrence.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* On the standardized 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, CURE records *Iris* $"Sil" = 0.489, "ARI" = 0.560$; *Customers* $"Sil" = 0.335$; *Multi-Shape* $"Sil" = 0.609$; and *High-D* $"Sil" = 0.238$, with mean runtime *$2.423$ ms*. Every number aligns tightly with our formal analysis:
- *Strong Convex & Elongated Recovery:* CURE achieves the highest *Customers* Silhouette ($0.335$) and ties CHAMELEON for the top *Iris* quality ($"Sil" = 0.489, "ARI" = 0.560$), outperforming BIRCH ($0.435 / 0.538$) by avoiding coarse spherical leaf quantization.
- *Multi-Shape Leap with Concavity Ceiling:* On *Multi-Shape*, CURE jumps to $"Sil" = 0.609$—substantially above BIRCH ($0.452$) and ROCK ($0.552$)—validating Theorem 2's anti-chaining multi-representative geometry, yet trailing CHAMELEON ($0.712$) due to the concentric/concave shrinkage limit proven in Proposition 1.
- *High-Dimensional Sensitivity & Speed:* At $2.423$ ms, CURE is $1.62 times$ faster than CHAMELEON ($3.922$ ms) and $10.5 times$ faster than ROCK ($25.428$ ms), while its *High-D* drop to $"Sil" = 0.238$ reflects Euclidean distance concentration weakening farthest-first boundary selection.
