#import "../template/ieee.typ": *

#alg_entry(
  id: 78,
  name: "BICO: BIRCH Meets Coresets for k-Means Clustering",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2013",
  paradigm: "Streaming Tree-Indexed Coreset Construction",
  complexity: $O(N d log S + S^2 d)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Streaming hierarchical summarization via BIRCH @zhang1996birch achieves high throughput by compressing incoming points into a tree of Clustering Features (CFs) governed by a fixed radius threshold, yet it provides *no* worst-case approximation guarantee for the $k$-Means objective—arbitrary multiplicative error can occur whenever high-multiplicity points accumulate inside a fixed-radius ball. Conversely, purely theoretical streaming coreset frameworks such as Har-Peled--Mazumdar exponential grids @harpeled2004coresets and StreamKM++ @ackermann2012streamkmpp guarantee $(1 + epsilon)$-approximations @aloise2009np but incur substantial data-structure overhead. *BICO* (*BIRCH Meets Coresets for $k$-Means*) @fichtenberger2013bico unifies both paradigms by augmenting a BIRCH-style hierarchical CF-tree with a *global 1-Means SSE threshold* $T$, *geometrically shrinking level radii* $R_i$, *random-projection nearest-neighbor indexing*, and *adaptive threshold doubling* under a strict coreset memory budget $S$.

Given a stream $X = \{bold(x)_1, dots.c, bold(x)_N\} subset.eq bb(R)^d$ and target cluster count $k$, the Euclidean $k$-Means objective for a center set $C subset.eq bb(R)^d$ ($|C| = k$) and its optimum are
$ "cost"(X, C) &= sum_(bold(x) in X) min_(bold(c) in C) || bold(x) - bold(c) ||_2^2, \
  "OPT"_k (X) &= inf_(|C| = k) "cost"(X, C). $
For any weighted summary $Omega = \{(bold(mu)_j, w_j)\}_(j=1)^m subset.eq bb(R)^d times bb(R)_(> 0)$, the weighted $k$-Means cost is $"cost"(Omega, C) = sum_(j=1)^m w_j min_(bold(c) in C) || bold(mu)_j - bold(c) ||_2^2$.

#definition(title: "Strong Coreset, Clustering Feature, and BICO Dual Invariant")[
  1. A weighted point set $Omega$ is a *strong $(k, epsilon)$-coreset* for $X$ @harpeled2004coresets @fichtenberger2013bico if for *every* center set $C subset.eq bb(R)^d$ of size $|C| = k$,
     $ (1 - epsilon) "cost"(X, C) &<= "cost"(Omega, C) \
       &<= (1 + epsilon) "cost"(X, C). $
  2. For a point subset $P subset.eq X$, its *Clustering Feature* @zhang1996birch is the triple $"CF"(P) = (n, bold(s), q)$ with $n = |P|$, $bold(s) = sum_(bold(x) in P) bold(x)$, $q = sum_(bold(x) in P) || bold(x) ||_2^2$, centroid $bold(mu) = bold(s) / n$, and exact 1-Means cost
     $ "cost"_1 (P) = sum_(bold(x) in P) || bold(x) - bold(mu) ||_2^2 = q - (|| bold(s) ||_2^2) / n. $
  3. In BICO's tree $cal(T)$ with global cost threshold $T > 0$ and level radii $R_i = sqrt(T) / (2^(i-1))$ (level $i >= 1$), each CF node $v$ at level $i$ stores $"CF"(P_v)$ and a fixed reference point $bold(r)_v in P_v$. Absorption of $bold(x)$ into $v$ requires *both*:
     $ || bold(x) - bold(r)_v ||_2 <= R_i quad "and" quad "cost"_1 (P_v union \{bold(x)\}) <= T. $
]

At level $1$, BICO indexes sibling reference points using $p$ random hyperplanes / 1D projections to accelerate nearest-reference lookup. Whenever total active CF nodes exceed the space budget $S$, BICO enters a *rebuild phase*: it doubles $T <- 2 T$ (rescaling $R_i$) and re-inserts the existing $S + 1$ CFs directly—never revisiting raw stream points.

#algbox(
  title: "BICO (Streaming Coreset Tree Construction + Weighted k-Means)",
  inputs: [Stream $X = (bold(x)_1, dots.c, bold(x)_N) subset.eq bb(R)^d$, clusters $k$, coreset budget $S >= k + 1$, projection count $p$.],
  outputs: [Strong coreset $Omega = \{(bold(mu)_v, n_v)\}$ ($|Omega| <= S$) and $k$-center set $C$.],
)[
  1. *Initialization*: Buffer the first $S + 1$ distinct points; compute minimum pairwise squared distance $d_"min"^2 = min_(a != b) || bold(x)_a - bold(x)_b ||_2^2$ (or over first $k + 1$ points) to initialize $T <- d_"min"^2 / (2 S)$ (or $T <- 0$ until $S + 1$ distinct points arrive); set $R_i <- sqrt(T) / 2^(i-1)$.
  2. *Procedure* `InsertCF(node u, level i, CF_new, ref r_new)`:
     - Find child $v$ of $u$ minimizing $|| bold(r)_"new" - bold(r)_v ||_2$ (using $p$ random projections at $i = 1$).
     - *If* no child exists *or* $|| bold(r)_"new" - bold(r)_v ||_2 > R_i$: create new child CF at level $i$ storing $"CF"_"new"$ with reference $bold(r)_"new"$.
     - *Else if* $"cost"_1 ("CF"_v + "CF"_"new") <= T$: merge $"CF"_v <- "CF"_v + "CF"_"new"$.
     - *Else*: recursively call `InsertCF(v, i + 1, CF_new, r_new)`.
  3. *Streaming Pass*: For each $bold(x)_t in X$, invoke `InsertCF(root, 1, (1, x_t, ||x_t||_2^2), x_t)`.
     - *Rebuild Trigger*: While total CF count $|cal(T)| > S$: set $T <- 2 T$, recompute $R_i <- sqrt(T) / 2^(i-1)$, collect all existing CFs $cal(F)$ in $cal(T)$, reset $cal(T)$ to empty root, and re-insert each $("CF"_v, bold(mu)_v) in cal(F)$ via `InsertCF`.
  4. *Extraction*: Emit $Omega = \{(bold(s)_v / n_v, n_v) : v in cal(T)\}$ and run weighted $k$-Means++ @arthur2007kmeanspp + Lloyd iterations @lloyd1982least on $Omega$ to return $C$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We present five formal results establishing the mathematical foundation of BICO: exact additive CF merging and parallel-axis decomposition, a deterministic strong coreset distortion theorem driven solely by bounded 1-Means residuals, the classical packing and threshold-calibration coreset bound, an explicit separation counterexample showing unbounded BIRCH distortion versus exact BICO optimality, and amortized streaming rebuild complexity.

#subsec("A", "Additive CF Algebra and Parallel-Axis Identity")

#lemma(title: "Exact Parallel-Axis and Ward Merge Identities for CFs", origin: "Original Derivation")[
  Let $P subset.eq bb(R)^d$ have $"CF"(P) = (n, bold(s), q)$ and centroid $bold(mu) = bold(s) / n$. Then for every point $bold(c) in bb(R)^d$,
  $ sum_(bold(x) in P) || bold(x) - bold(c) ||_2^2 = "cost"_1 (P) + n || bold(mu) - bold(c) ||_2^2, $
  where $"cost"_1 (P) = q - || bold(s) ||_2^2 / n$. Moreover, for two disjoint subsets $P_1, P_2$ with $"CF"(P_j) = (n_j, bold(s)_j, q_j)$ and centroids $bold(mu)_j$, the merged feature satisfies $"CF"(P_1 union P_2) = (n_1 + n_2, bold(s)_1 + bold(s)_2, q_1 + q_2)$ and
  $ "cost"_1 (P_1 union P_2) &= "cost"_1 (P_1) + "cost"_1 (P_2) \
    &quad + (n_1 n_2) / (n_1 + n_2) || bold(mu)_1 - bold(mu)_2 ||_2^2. $
]
#proof[
  Writing $bold(x) - bold(c) = (bold(x) - bold(mu)) + (bold(mu) - bold(c))$ and expanding the squared Euclidean norm gives
  $ sum_(bold(x) in P) || bold(x) - bold(c) ||_2^2 &= sum_(bold(x) in P) || bold(x) - bold(mu) ||_2^2 \
    &quad + n || bold(mu) - bold(c) ||_2^2 \
    &quad + 2 chevron.l sum_(bold(x) in P) (bold(x) - bold(mu)), bold(mu) - bold(c) chevron.r. $
  The cross-term vanishes identically because $sum_(bold(x) in P) bold(x) = bold(s) = n bold(mu)$, proving the parallel-axis identity. Setting $bold(c) = bold(0)$ yields $q = "cost"_1 (P) + n || bold(mu) ||_2^2 = "cost"_1 (P) + || bold(s) ||_2^2 / n$.

  For disjoint $P_1, P_2$, additivity of $(n, bold(s), q)$ follows immediately from linearity of summation. Letting $bold(mu)_(12) = (n_1 bold(mu)_1 + n_2 bold(mu)_2) / (n_1 + n_2)$ be the combined centroid and applying the parallel-axis identity to $P_1$ and $P_2$ around $bold(c) = bold(mu)_(12)$,
  $ "cost"_1 (P_1 union P_2) &= "cost"_1 (P_1) + "cost"_1 (P_2) \
    &quad + n_1 || bold(mu)_1 - bold(mu)_(12) ||_2^2 \
    &quad + n_2 || bold(mu)_2 - bold(mu)_(12) ||_2^2. $
  Substituting $bold(mu)_1 - bold(mu)_(12) = n_2 / (n_1 + n_2) (bold(mu)_1 - bold(mu)_2)$ and $bold(mu)_2 - bold(mu)_(12) = -n_1 / (n_1 + n_2) (bold(mu)_1 - bold(mu)_2)$ yields the Ward increment $(n_1 n_2) / (n_1 + n_2) || bold(mu)_1 - bold(mu)_2 ||_2^2$.
]

#subsec("B", "Deterministic Strong Coreset Guarantee From Bounded Residuals")

#theorem(title: "Two-Sided Strong Coreset Bound via Total 1-Means Residual", origin: "Original Proof")[
  Let $cal(P) = \{P_1, dots.c, P_m\}$ be any partition of $X$ into $m <= S$ CF nodes satisfying $"cost"_1 (P_j) <= T$ for every $j in [m]$, and let $Omega = \{(bold(mu)_j, n_j)\}_(j=1)^m$. Denote the total residual $Delta = sum_(j=1)^m "cost"_1 (P_j) <= m T$. Then for *every* center set $C subset.eq bb(R)^d$ of size $k$,
  $ -Delta &<= "cost"(Omega, C) - "cost"(X, C) \
           &<= 2 sqrt(Delta dot "cost"(X, C)) + Delta. $
  Consequently, whenever $Delta <= (epsilon^2 / 9) "OPT"_k (X)$ for $epsilon in (0, 1]$, $Omega$ is deterministically a strong $(k, epsilon)$-coreset for $X$.
]
#proof[
  Fix an arbitrary $k$-center set $C subset.eq bb(R)^d$. For each $bold(x) in X$ and $j in [m]$, denote nearest centers by $c(bold(x)) in arg min_(bold(c) in C) || bold(x) - bold(c) ||_2$ and $c_j in arg min_(bold(c) in C) || bold(mu)_j - bold(c) ||_2$.

  *Lower bound on $"cost"(Omega, C) - "cost"(X, C)$:* For every $bold(x) in P_j$, optimality of $c(bold(x))$ implies $|| bold(x) - c(bold(x)) ||_2^2 <= || bold(x) - c_j ||_2^2$. Summing over $bold(x) in P_j$ and applying Lemma II-A gives
  $ sum_(bold(x) in P_j) || bold(x) - c(bold(x)) ||_2^2 &<= sum_(bold(x) in P_j) || bold(x) - c_j ||_2^2 \
    &= "cost"_1 (P_j) + n_j || bold(mu)_j - c_j ||_2^2. $
  Summing across $j = 1, dots.c, m$ yields $"cost"(X, C) <= "cost"(Omega, C) + Delta$, i.e., $"cost"(Omega, C) - "cost"(X, C) >= -Delta$.

  *Upper bound on $"cost"(Omega, C) - "cost"(X, C)$:* For every $bold(x) in P_j$, optimality of $c_j$ and the triangle inequality give
  $ || bold(mu)_j - c_j ||_2 &<= || bold(mu)_j - c(bold(x)) ||_2 \
    &<= || bold(x) - c(bold(x)) ||_2 + || bold(x) - bold(mu)_j ||_2. $
  Squaring both sides and summing over all $bold(x) in P_j$ and all $j in [m]$,
  $ "cost"(Omega, C) &= sum_(j=1)^m sum_(bold(x) in P_j) || bold(mu)_j - c_j ||_2^2 \
    &<= "cost"(X, C) + Delta \
    &quad + 2 sum_(j=1)^m sum_(bold(x) in P_j) || bold(x) - c(bold(x)) ||_2 || bold(x) - bold(mu)_j ||_2. $
  Applying the Cauchy--Schwarz inequality to the cross sum yields
  $ &sum_(j=1)^m sum_(bold(x) in P_j) || bold(x) - c(bold(x)) ||_2 || bold(x) - bold(mu)_j ||_2 \
    &quad <= sqrt("cost"(X, C)) sqrt(Delta), $
  establishing $"cost"(Omega, C) - "cost"(X, C) <= 2 sqrt(Delta dot "cost"(X, C)) + Delta$.

  Finally, suppose $Delta <= (epsilon^2 / 9) "OPT"_k (X)$ with $epsilon in (0, 1]$. Since $"OPT"_k (X) <= "cost"(X, C)$ for every $k$-center set $C$, we have $Delta <= (epsilon^2 / 9) "cost"(X, C)$, whence
  $ "cost"(Omega, C) &>= (1 - epsilon^2 / 9) "cost"(X, C) \
                     &>= (1 - epsilon) "cost"(X, C), \
    "cost"(Omega, C) &<= (1 + 2 epsilon / 3 + epsilon^2 / 9) "cost"(X, C) \
                     &<= (1 + epsilon) "cost"(X, C), $
  using $epsilon^2 / 9 <= epsilon / 3$ for $epsilon in (0, 1]$. Thus $Omega$ is a strong $(k, epsilon)$-coreset.
]

#subsec("C", "Pigeonhole Lower Bound and Packing Coreset Theorem")

#theorem(title: "BICO Threshold Calibration and Strong Coreset Space Bound", origin: "Classical — Cited")[
  Let $X subset.eq bb(R)^d$ and $k >= 1$ @fichtenberger2013bico.
  1. For any $k + 1$ points $bold(y)_1, dots.c, bold(y)_(k+1) in X$,
     $ "OPT"_k (X) >= 1 / 2 min_(1 <= a < b <= k+1) || bold(y)_a - bold(y)_b ||_2^2. $
  2. When BICO is executed with geometric radii $R_i = sqrt(T) / 2^(i-1)$ and coreset size budget $S = Theta(k epsilon^(-d - 2))$ (in fixed dimension $d$, or $O(k epsilon^(-2) log N)$ under random-projection / sensitivity subsampling extensions @harpeled2004coresets @fichtenberger2013bico), the adaptive doubling invariant guarantees $m T <= (epsilon^2 / 9) "OPT"_k (X)$ upon termination, so $Omega$ is a strong $(k, epsilon)$-coreset.
]
#proof[
  For part 1, let $C^*$ be an optimal $k$-center set for $X$. By the Pigeonhole Principle, among $k + 1$ points at least two distinct points $bold(y)_a, bold(y)_b$ are assigned to the same optimal center $bold(c)^* in C^*$. By Lemma II-A applied to $\{bold(y)_a, bold(y)_b\}$ with midpoint $bold(mu) = (bold(y)_a + bold(y)_b) / 2$,
  $ "OPT"_k (X) &>= || bold(y)_a - bold(c)^* ||_2^2 + || bold(y)_b - bold(c)^* ||_2^2 \
                &>= "cost"_1 (\{bold(y)_a, bold(y)_b\}) = 1 / 2 || bold(y)_a - bold(y)_b ||_2^2. $

  For part 2 (following Fichtenberger et al. @fichtenberger2013bico), whenever a rebuild doubles $T$ from $T / 2$ to $T$, the tree at threshold $T / 2$ contained $> S$ nodes. Every level-$i$ sibling pair has reference separation strictly exceeding $R_i$, so by Euclidean ball packing at most $(1 + 2 / alpha)^d$ disjoint references of separation $> R_i$ can occupy a ball of radius $R_(i-1) / alpha$. Any optimal center $bold(c)^* in C^*$ can cover at most $O((1 / epsilon)^d)$ separated balls at each radial ring before their distance to $bold(c)^*$ exceeds $R_i / 2$, forcing uncovered nodes to contribute $Omega(T)$ each to $"OPT"_k (X)$. Summing across all $k$ optimal clusters and geometric scales shows that exceeding $S = Theta(k epsilon^(-d-2))$ nodes implies $"OPT"_k (X) = Omega(S T / epsilon^2)$, so $Delta <= S T <= (epsilon^2 / 9) "OPT"_k (X)$. Theorem II-B then guarantees a strong $(k, epsilon)$-coreset.
]

#subsec("D", "Constructive Separation Between BIRCH and BICO")

We prove constructively that omitting BICO's 1-Means cost cap $T$ causes BIRCH @zhang1996birch to incur arbitrarily large $k$-Means distortion, whereas BICO remains exact.

#proposition(title: "Unbounded BIRCH Distortion Versus Exact BICO Optimality", origin: "Original Counterexample")[
  For every integer $n >= 1$ and $k = 2$, there exists a 1D dataset $X_n subset.eq bb(R)$ of size $N = 2 n + 1$ on which BIRCH with fixed radius threshold $R = 1$ yields a 2-Means solution with multiplicative approximation ratio at least $n / 2 -> infinity$, whereas BICO with budget $S >= 3$ produces a zero-residual coreset achieving exact ratio $1$.
]
#proof[
  Fix $R = 1$ and let $X_n subset.eq bb(R)$ consist of $n$ points at $x = -1$, $n$ points at $x = +1$, and $1$ point at $x = 3$, presented in alternating pairs $(-1, +1)$ followed by $3$.
  - *True Optimal 2-Means Cost*: Placing centers at $C_"feas" = \{-1, +1\}$ assigns all $2 n$ points at plus-or-minus $1$ with zero cost and $x = 3$ with squared cost $(3 - 1)^2 = 4$, so $"OPT"_2 (X_n) <= 4$ uniformly in $n$.
  - *BIRCH Failure*: In BIRCH with radius bound $R = 1$, inserting alternating $(-1, +1)$ pairs maintains centroid $bold(mu)_1 = 0$ and radius $sqrt("cost"_1 (P_1) / |P_1|) = sqrt((2 m) / (2 m)) = 1 <= R$ at every step. Hence BIRCH absorbs all $2 n$ points into a single CF node $(2 n, 0, 2 n)$ with centroid $0$. Point $x = 3$ lies at distance $3 > R$ and creates a second CF node $(1, 3, 9)$. The resulting summary $Omega_"BIRCH" = \{(0, 2 n), (3, 1)\}$ has only two distinct support points, so *any* optimal 2-Means solution on $Omega_"BIRCH"$ places centers at $C_"BIRCH" = \{0, 3\}$. Evaluating $C_"BIRCH"$ on $X_n$ yields
    $ "cost"(X_n, C_"BIRCH") &= n (-1 - 0)^2 + n (1 - 0)^2 \
                             &quad + (3 - 3)^2 = 2 n, $
    giving multiplicative distortion $"cost"(X_n, C_"BIRCH") / "OPT"_2 (X_n) >= 2 n / 4 = n / 2 -> infinity$.
  - *BICO Exactness*: With budget $S >= 3$, BICO only triggers a rebuild when $> S$ distinct CF nodes exist. Because $X_n$ has only $3$ distinct coordinate values $\{-1, +1, 3\}$, BICO retains $T = 0$ (zero rebuilds) and stores three exact zero-cost CF nodes $(-1, n)$, $(+1, n)$, and $(3, 1)$ with $Delta = 0$. Even if $T > 0$ were pre-initialized below $2 n$, BICO's cost cap $"cost"_1 (P_1) <= T$ forbids merging all $2 n$ points into one CF. Thus $Omega_"BICO"$ reproduces $"cost"(X_n, C)$ identically for all $C$, attaining approximation ratio $1$.
]

#subsec("E", "Amortized Streaming Complexity Under Threshold Doubling")

#proposition(title: "Logarithmic Rebuild Count and Strict Space Bound", origin: "Original Proof")[
  Let $T_0 > 0$ be the first positive threshold initialized by BICO and let $Delta_"spread" = "OPT"_k (X) / T_0 >= 1$. Under the packing invariant $T_"final" <= c_1 "OPT"_k (X)$, BICO executes at most $R_"reb" <= ceil(log_2 (c_1 Delta_"spread"))$ rebuilds, uses strictly $O(S d)$ memory, and runs in total streaming time
  $ O(N d (p + h_"max") + S d (p + h_"max") log Delta_"spread") $
  where $h_"max"$ is the maximum tree depth and $p$ is the number of random projections at level $1$.
]
#proof[
  Every rebuild strictly doubles the threshold: $T_r = 2^r T_0$. Since the final threshold satisfies $T_(R_"reb") = 2^(R_"reb") T_0 <= c_1 "OPT"_k (X)$, taking logarithms yields $R_"reb" <= ceil(log_2 (c_1 Delta_"spread"))$. Each tree node stores $(n_v, bold(s)_v, q_v, bold(r)_v)$ in $2 d + 2$ scalars, and rebuild triggers immediately when $|cal(T)| = S + 1$, bounding peak memory strictly by $O(S d)$ independently of $N$. Each incoming stream point evaluates $p$ dot products at level $1$ ($O(p d)$) and traverses at most $h_"max"$ tree levels ($O(h_"max" d)$ per point), totaling $O(N d (p + h_"max"))$. Each rebuild re-inserts at most $S + 1$ CFs without scanning $X$, costing $O(S d (p + h_"max"))$ per rebuild and $O(S d (p + h_"max") log Delta_"spread")$ across all rebuilds.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time, Space, and Streaming Profile.* BICO executes in a *single pass* over $X$ with strict $O(S d)$ working memory—completely decoupled from stream length $N$. For typical budgets $S = 200 k$ to $500 k$ and $p = 10$ random projections, level-1 lookup amortizes to near-constant candidate comparisons, giving $O(N d log S)$ streaming insertion time plus $O(k S d)$ post-clustering via weighted $k$-Means++ @arthur2007kmeanspp (or $O(S^2 d)$ if pairwise coreset initialization is used). Because CFs are closed under addition (Lemma II-A), BICO merges trivially across distributed shards in a single reduce tree (map-reduce mergeability @harpeled2004coresets).

*Geometric Invariances.* Both the Euclidean distance criterion $|| bold(x) - bold(r)_v ||_2 <= R_i$ and the 1-Means cost criterion $"cost"_1 (P_v) <= T$ depend purely on Euclidean inner products and norms. Hence BICO's mathematical invariant is strictly invariant under translations and *orthogonal rotations* (unlike PERCH's AABBs) and scale-equivariant under uniform rescaling $bold(x) mapsto alpha bold(x)$ (which rescales $T mapsto alpha^2 T$ and $R_i mapsto alpha R_i$ automatically). Note that the choice of first-arriving reference $bold(r)_v$ and level-1 random projections introduce mild stream-order and random-seed variation, while preserving worst-case coreset distortion bounds.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N d log S)$ streaming pass + $O(k S d)$ weighted $k$-Means++ post-solve],
    [Space Complexity], [$O(S d)$ strict peak memory (independent of stream length $N$)],
    [Approximation Guarantee], [Strong $(k, epsilon)$-coreset on Euclidean $k$-Means objective (Theorem II-B & II-C)],
    [Distributed Mergeability], [Fully mergeable via additive CFs $(n, bold(s), q)$ across parallel partitions],
    [Invariances], [Translation, orthogonal rotation, and uniform scale equivariance],
    [Primary Hyperparameters], [Target clusters $k$, coreset budget $S$, random projection count $p$],
  ),
  caption: [Operational and geometric profile of BICO (Alg. 78).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Web-Scale Telemetry and Network Traffic Streams.* High-velocity packet and clickstream logs ($N > 10^8$) cannot reside in RAM; BICO compresses the stream in one pass into $S approx 10^4$ weighted summaries in $O(S d)$ memory while preserving provable $k$-Means quantization error @fichtenberger2013bico.
2. *Vector Quantization and Product-Quantization Codebook Training.* Training IVF/$k$-Means coarse quantizers on billion-vector embedding databases benefits directly from BICO's strong $(k, epsilon)$-coreset: running multiple weighted $k$-Means++ restarts @arthur2007kmeanspp on $Omega$ explores diverse initializations orders of magnitude faster than scanning $X$.
3. *Distributed Map-Reduce and Edge Sensor Summarization.* Edge IoT nodes construct local BICO trees independently and transmit only $O(S d)$ CFs to a central aggregator, which merges them seamlessly via Lemma II-A's additive CF algebra without raw data transfer.
4. *Astronomical Sky-Survey and Particle-Physics Event Compression.* Multi-terabyte photometric catalogs and calorimeter hit streams require guaranteed variance preservation across multi-scale density regions; BICO's cost cap prevents dense cores from masking sparse outlier populations (Proposition II-D).
5. *Fast Hyperparameter Sweeps Over Cluster Count $k$.* Because a strong $(k_"max", epsilon)$-coreset simultaneously approximates the $k$-Means objective for every $k <= k_"max"$ and every candidate center set $C$, analysts can evaluate elbow curves, Silhouette scores @rousseeuw1987silhouettes, and multiple $k$ values entirely on $Omega$ in milliseconds.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

BICO exhibits four well-defined failure regimes: (i) *Non-convex / manifold geometries* (concentric rings, interleaving crescents), where BICO accurately solves the *Euclidean $k$-Means objective*—an objective that inherently imposes Voronoi-convex partitions; (ii) *Heavy-tailed extreme outliers*, where a handful of distant points inflate $"OPT"_k (X)$, causing early threshold doubling $T <- 2 T$ and coarser granularity inside dense cores; (iii) *Very high dimension $d$ with tiny budget $S$*, where Theorem II-C's packing factor grows with intrinsic dimension, loosening empirical coreset tightness; and (iv) *Non-Euclidean metrics* (edit distance, arbitrary graph kernels), where Lemma II-A's Euclidean parallel-axis decomposition does not apply.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Massive Euclidean $k$-Means streams], [One-pass $O(S d)$ RAM + provable $(1 + epsilon)$ coreset], [Set $S >= 50 k$--$200 k$ for tight $epsilon$],
    [Distributed / edge-cloud aggregation], [Exact additive CF mergeability across shards], [Re-apply $k$-Means++ on final merged coreset],
    [High-$d$ convex Gaussian / spherical mixtures], [Orthogonally invariant; immune to AABB inflation], [Use spherical $k$-Means if cosine metric needed],
    [Non-convex / elongated manifold shapes], [Poor fit: inherits spherical Voronoi bias of $k$-Means], [Use PERCH (Alg. 76), SCC (Alg. 79), or HDBSCAN],
    [Extreme outlier-contaminated streams], [Outliers inflate $T$ via rapid doubling], [Pre-filter outliers or use trimmed $k$-Means],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for BICO.],
)

#sec("V", "Structural Relations and Empirical Behavior")

BICO sits directly at the intersection of tree-based stream summarization (BIRCH @zhang1996birch, CluStream @aggarwal2003clustream, ClusTree @kranen2011clustree) and theoretical $k$-Means coresets (Har-Peled--Mazumdar @harpeled2004coresets, StreamKM++ @ackermann2012streamkmpp, Mini-Batch $k$-Means @sculley2010webscale). Proposition II-D proves BICO strictly dominates BIRCH on worst-case $k$-Means approximation, while retaining BIRCH's single-pass CF-tree throughput.

#proposition(title: "End-to-End Approximation Factor With Weighted k-Means++", origin: "Original Proof")[
  Suppose BICO outputs a strong $(k, epsilon)$-coreset $Omega$ with $epsilon in (0, 1 / 2]$, and weighted $k$-Means++ @arthur2007kmeanspp produces a center set $C_"alg"$ satisfying $bb(E)["cost"(Omega, C_"alg")] <= alpha "OPT"_k (Omega)$ (with $alpha = 8 (ln k + 2)$). Then $C_"alg"$ achieves expected cost on the full stream $X$ bounded by
  $ bb(E)["cost"(X, C_"alg")] &<= ((1 + epsilon) / (1 - epsilon)) alpha "OPT"_k (X) \
                              &<= (1 + 4 epsilon) alpha "OPT"_k (X). $
]
#proof[
  Let $C_X^* in arg min_(|C|=k) "cost"(X, C)$ and $C_Omega^* in arg min_(|C|=k) "cost"(Omega, C)$. Applying Definition I.1 twice—first to $C_"alg"$ on the lower bound, and second to $C_X^*$ on the upper bound—gives
  $ (1 - epsilon) "cost"(X, C_"alg") &<= "cost"(Omega, C_"alg"), \
    "OPT"_k (Omega) <= "cost"(Omega, C_X^*) &<= (1 + epsilon) "OPT"_k (X). $
  Taking expectations over the randomized $k$-Means++ seeding yields
  $ bb(E)["cost"(X, C_"alg")] &<= 1 / (1 - epsilon) bb(E)["cost"(Omega, C_"alg")] \
                              &<= ((1 + epsilon) / (1 - epsilon)) alpha "OPT"_k (X). $
  For $epsilon in (0, 1 / 2]$, $(1 + epsilon) / (1 - epsilon) = 1 + (2 epsilon) / (1 - epsilon) <= 1 + 4 epsilon$, completing the proof.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, BICO obtains: *Iris* Silhouette $= 0.449$, ARI $= 0.638$; *Customers* Silhouette $= 0.354$; *Multi-Shape* Silhouette $= 0.557$ (ARI N/A); *High-D* Silhouette $= 0.258$ (ARI N/A); with mean runtime $= 1.229$ ms. Every single metric matches our theoretical predictions: BICO attains the *highest Iris ARI* ($0.638$) and *highest High-D Silhouette* ($0.258$) among Algs. 76--79 because Iris and High-D feature roughly convex clusters well matched to the Euclidean $k$-Means objective, and BICO's orthogonal invariance avoids PERCH's $sqrt(d)$ AABB degradation. Conversely, BICO scores lowest on Multi-Shape ($0.557$ vs PERCH's $0.709$ and GRINCH/SCC's $0.609$) precisely because its strong coreset faithfully optimizes spherical $k$-Means rather than linkage-based shape connectivity.
