#import "../template/ieee.typ": *

#alg_entry(
  id: 9,
  name: "Bisecting K-Means",
  category: "Cat. I: Centroid & Partitioning",
  year: "2000",
  paradigm: "Divisive Hierarchical Binary Prototype Splitting",
  complexity: $O(n d R T log k)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Flat $k$-means clustering (Alg. 1) optimizes all $k$ prototypes simultaneously and therefore suffers from severe sensitivity to $k$-way initialization, frequent empty clusters in high dimension, and linear $O(k)$ distance evaluations per point per iteration. *Bisecting K-Means*, introduced and systematically benchmarked by Steinbach, Karypis, and Kumar @steinbach2000comparison for large sparse document corpora, replaces simultaneous $k$-way partitioning with a top-down *divisive binary hierarchy* built solely from repeated 2-means splits.

Given a finite dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$ and target leaf count $k in {1, dots, n}$, define for any non-empty subset $C subset.eq {1, dots, n}$ its arithmetic centroid $overline(x)_C = |C|^(-1) sum_(i in C) x_i$ and within-cluster sum of squares
$ "SSE"(C) = sum_(i in C) norm(x_i - overline(x)_C)_2^2. $
Starting from the root partition $cal(C)^((1)) = { {1, dots, n} }$, Bisecting K-Means performs $k - 1$ sequential binary splits. At step $m in {1, dots, k - 1}$, given current leaf clusters $cal(C)^((m)) = {C_1, dots, C_m}$, the algorithm: (i) selects one active cluster $C_(j^*)$ according to a priority rule—canonically maximum $"SSE"(C_j)$ (or alternatively largest cardinality $|C_j|$); (ii) solves $R >= 1$ independent trials of Euclidean 2-means on $C_(j^*)$ and retains the best bisection $(A^*, B^*)$ minimizing $"SSE"(A) + "SSE"(B)$; and (iii) replaces $C_(j^*)$ with $A^*, B^*$ to form $cal(C)^((m+1))$.

#definition(title: "Divisive Binary Centroid Hierarchy")[
  A *bisecting $k$-means hierarchy* on $cal(X)$ is a rooted full binary tree $cal(T)$ with $k$ leaves and $k - 1$ internal nodes such that: (i) the root represents ${1, dots, n}$; (ii) each internal node $u$ stores cluster $C_u$, centroid $overline(x)_(C_u)$, and child pair $(C_(u_L), C_(u_R))$ constituting a discrete 2-means CVT separated by the perpendicular bisector hyperplane
  $ H_u = { x in bb(R)^d : &chevron.l x - 1 / 2 (mu_(u_L) + mu_(u_R)), \
    &space mu_(u_L) - mu_(u_R) chevron.r = 0 }; $
  and (iii) the leaf family $cal(C)_"leaf" (cal(T))$ forms a disjoint $k$-partition of ${1, dots, n}$ with total inertia $J^*(cal(C)_"leaf") = sum_(v in "leaves"(cal(T))) "SSE"(C_v)$.
]

Algorithm 1 formalizes the canonical max-SSE multi-trial procedure together with optional global Lloyd post-refinement.

#algbox(
  title: "Bisecting K-Means Algorithm",
  inputs: [Dataset $cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$, target clusters $k >= 1$, trials per split $R >= 1$, 2-means cap $T$.],
  outputs: [Partition $cal(C) = {C_1, dots, C_k}$, centroids $bold(mu)$, binary tree $cal(T)$.],
)[
  1. Initialize leaf pool $cal(P) = { {1, dots, n} }$ in a max-heap keyed by $"SSE"(C)$, and initialize tree $cal(T)$ with root node $r$.
  2. *While* $|cal(P)| < k$ *do*:
  3. #h(1em) Pop cluster $C$ maximizing $"SSE"(C)$ (if $"SSE"(C) = 0$, *break*).
  4. #h(1em) *For* trial $tau = 1, dots, R$ *do*:
  5. #h(2em) Seed two distinct centroids $(mu_1^((0)), mu_2^((0)))$ from $C$ (e.g., via K-Means++ @arthur2007kmeanspp) and run Lloyd's 2-means (Alg. 1) on $C$ for at most $T$ steps to obtain $(A_tau, B_tau)$.
  6. #h(1em) Select $(A^*, B^*) = arg min_(1 <= tau <= R) [ "SSE"(A_tau) + "SSE"(B_tau) ]$.
  7. #h(1em) Attach $(A^*, B^*)$ as children of $C$ in $cal(T)$ and push $A^*, B^*$ onto $cal(P)$.
  8. *Optional refinement:* Run global Lloyd K-Means from the $k$ leaf centroids until CVT convergence.
  9. *Return* leaf partition $cal(C)$, centroids $bold(mu)$, and hierarchy $cal(T)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five foundational results governing Bisecting K-Means: (1) an exact variance-reduction identity and quantitative lower bound proving strict energy decrease per split; (2) an explicit counterexample proving that even with an *exact global 2-means oracle*, Bisecting K-Means suffers an unbounded $k$-means approximation ratio and produces non-Voronoi partitions; (3) exact tree-height runtime bounds contrasting $O(n d R T log k)$ balanced performance against $O(n k d R T)$ worst-case degradation; (4) local 2-means NP-hardness and Dasgupta hierarchical cost behavior; and (5) geometric invariances and a subtle outlier-peeling breakdown characterization.

#subsec("A", "Exact Variance-Reduction Identity and Split-Gain Lower Bound")

#lemma(title: "Exact Binary Split Variance-Reduction Identity", origin: "Original Proof")[
  For any non-empty cluster $C subset.eq {1, dots, n}$ partitioned into disjoint non-empty subsets $A, B$ with centroids $mu_A, mu_B$, the reduction in within-cluster sum of squares satisfies identically
  $ Delta "SSE"(C -> A, B) &:= "SSE"(C) - "SSE"(A) - "SSE"(B) \
    &= (|A| |B|) / (|A| + |B|) norm(mu_A - mu_B)_2^2. $
]

#proof[
  Let $mu_C = (|A| mu_A + |B| mu_B) / |C|$ with $|C| = |A| + |B|$. Applying the parallel-axis identity separately to $A$ and $B$ relative to $mu_C$ yields
  $ "SSE"(C) &= sum_(i in A) norm(x_i - mu_C)_2^2 + sum_(i in B) norm(x_i - mu_C)_2^2 \
    &= "SSE"(A) + |A| norm(mu_A - mu_C)_2^2 + "SSE"(B) + |B| norm(mu_B - mu_C)_2^2. $
  Substituting $mu_A - mu_C = (|B| / |C|)(mu_A - mu_B)$ and $mu_B - mu_C = -(|A| / |C|)(mu_A - mu_B)$ gives
  $ Delta "SSE"(C -> A, B) &= (|A| |B|^2 + |B| |A|^2) / |C|^2 norm(mu_A - mu_B)_2^2 \
    &= (|A| |B|) / (|A| + |B|) norm(mu_A - mu_B)_2^2, $
  completing the identity.
]

#theorem(title: "Strict Monotone Descent and Quantitative Split-Gain Guarantee", origin: "Original Proof")[
  At any step $m in {1, dots, k - 1}$ of Algorithm 1, suppose cluster $C$ contains at least two distinct points and at least one trial seeds $(mu_1^((0)), mu_2^((0))) = (x_p, x_q)$ with $x_p != x_q in C$. Then the converged bisection $(A^*, B^*)$ produces two non-empty clusters and strictly lowers total inertia by at least
  $ J^*(cal(C)^((m))) - J^*(cal(C)^((m+1))) >= 1 / (|C|) norm(x_p - x_q)_2^2 > 0. $
  Furthermore, if a trial seeds $(overline(x)_C, x_(i^*))$ where $x_(i^*) in arg max_(i in C) norm(x_i - overline(x)_C)_2^2$, then $Delta "SSE"(C -> A^*, B^*) >= norm(x_(i^*) - overline(x)_C)_2^2 >= "SSE"(C) / |C|$.
]

#proof[
  Let $w = x_p - x_q != bold(0)$ and $tau = (norm(x_p)_2^2 - norm(x_q)_2^2) / 2$. In the initial Voronoi assignment $(A^((0)), B^((0)))$ of the trial seeded at $(x_p, x_q)$, every $i in A^((0))$ satisfies $norm(x_i - x_p)_2^2 <= norm(x_i - x_q)_2^2$, equivalently $chevron.l x_i, w chevron.r >= tau$, with strict margin $chevron.l x_p, w chevron.r - tau = norm(w)_2^2 / 2 > 0$ at $i = p$. Symmetrically, every $i in B^((0))$ satisfies $chevron.l x_i, w chevron.r <= tau$, with strict margin $chevron.l x_q, w chevron.r - tau = -norm(w)_2^2 / 2 < 0$ at $i = q$. Hence $p in A^((0))$ and $q in B^((0))$, so both sets are non-empty. Averaging over $A^((0))$ and $B^((0))$ yields
  $ chevron.l mu_(A^((0))), w chevron.r &>= tau + norm(w)_2^2 / (2 |A^((0))|), \
    chevron.l mu_(B^((0))), w chevron.r &<= tau - norm(w)_2^2 / (2 |B^((0))|). $
  Subtracting these two inequalities and applying Cauchy–Schwarz gives
  $ norm(mu_(A^((0))) - mu_(B^((0))))_2 norm(w)_2 &>= chevron.l mu_(A^((0))) - mu_(B^((0))), w chevron.r \
    &>= norm(w)_2^2 / 2 ( 1 / (|A^((0))|) + 1 / (|B^((0))|) ) \
    &= (|C| norm(w)_2^2) / (2 |A^((0))| |B^((0))|). $
  Dividing by $norm(w)_2 > 0$, substituting into Lemma 1, and invoking AM–GM ($|A^((0))| |B^((0))| <= |C|^2 / 4$) yields
  $ Delta "SSE"(C -> A^((0)), B^((0))) &= (|A^((0))| |B^((0))|) / |C| norm(mu_(A^((0))) - mu_(B^((0))))_2^2 \
    &>= |C| / (4 |A^((0))| |B^((0))|) norm(w)_2^2 >= 1 / (|C|) norm(x_p - x_q)_2^2. $
  Because subsequent Lloyd iterations on $C$ and multi-trial minimization over $tau = 1, dots, R$ never increase $"SSE"(A) + "SSE"(B)$, the accepted split $(A^*, B^*)$ satisfies the same lower bound.

  For the second bound, initializing at $(overline(x)_C, x_(i^*))$ produces initial assignment cost at most $sum_(i in C backslash {i^*}) norm(x_i - overline(x)_C)_2^2 + 0 = "SSE"(C) - norm(x_(i^*) - overline(x)_C)_2^2$; Lloyd monotonicity and the inequality $max_(i in C) norm(x_i - overline(x)_C)_2^2 >= |C|^(-1) "SSE"(C)$ complete the proof.
]

#subsec("B", "Greedy Irreversibility, Unbounded Ratio, and Non-Voronoi Leaves")

Unlike flat K-Means, Bisecting K-Means cannot undo an earlier hyperplane cut during subsequent splits. The next original counterexample proves that this hierarchical commitment causes arbitrarily large suboptimality—even when every local 2-means subproblem is solved to *exact global optimality*.

#theorem(title: "Unbounded Approximation Ratio and Non-Voronoi Leaves", origin: "Original Counterexample")[
  For any target factor $M > 0$ and $k = 3$, there exists a one-dimensional dataset $cal(X) subset.eq bb(R)$ on which Algorithm 1 without global post-refinement—even when every binary split solves Euclidean 2-means *globally to exact optimality*—produces a partition $cal(C)_"bi"$ satisfying:
  1. $J^*(cal(C)_"bi") / J^*(cal(C)^*) > M$, where $cal(C)^*$ is the global 3-means minimizer; and
  2. $cal(C)_"bi"$ is *not* a Voronoi equilibrium of its own leaf centroids.
]

#proof[
  Fix an integer $w >= max(3, ceil(sqrt(2 M)))$, set $L = w + 2$, and construct the symmetric 1D dataset $cal(X)$ of size $n = 2 w + 2$ consisting of $w$ points at $-L$, one point at $-1$, one point at $+1$, and $w$ points at $+L$ (perturbing co-located points by $epsilon -> 0^+$ yields distinct points while preserving all strict inequalities below).

  In one dimension, every optimal 2-means partition is a contiguous threshold cut separating distinct coordinates. By symmetry, only two candidate cuts exist at Step 1:
  - *Symmetric cut:* $A = {-L (times w), -1}$, $B = {+1, +L (times w)}$. Their centroids are $mu_A = -(w + 1)$ and $mu_B = (1 + w(w + 2)) / (w + 1) = w + 1$. By Lemma 1 on the atoms ${+1}$ and ${+L (times w)}$ separated by $L - 1 = w + 1$,
    $ "SSE"(A) = "SSE"(B) = (1 dot w) / (w + 1) (w + 1)^2 = w (w + 1), $
    giving total 2-means cost $"SSE"_"sym" = 2 w^2 + 2 w$.
  - *Asymmetric cut:* $A' = {-L (times w)}$, $B' = {-1, +1, +L (times w)}$. Here $"SSE"(A') = 0$ and $mu_(B') = (w (w + 2)) / (w + 2) = w$, so
    $ "SSE"_"asym" &= "SSE"(B') = (-1 - w)^2 + (1 - w)^2 + w (2)^2 \
      &= 2 w^2 + 4 w + 2. $
  Because $"SSE"_"sym" = 2 w^2 + 2 w < 2 w^2 + 4 w + 2 = "SSE"_"asym"$ for every $w >= 1$, the symmetric split $(A, B)$ is the *unique global minimizer* of 2-means on $cal(X)$.

  At Step 2 ($k = 3$), both $A$ and $B$ have identical size $w + 1$ and identical $"SSE" = w (w + 1)$. Hence under max-SSE, max-cardinality, or max-gain selection, Algorithm 1 splits one branch (without loss of generality $A$) into ${-L (times w)}$ and ${-1}$ (zero SSE each) while leaving $B$ intact, yielding
  $ cal(C)_"bi" = { {-L (times w)}, {-1}, {+1, +L (times w)} } $
  with total cost $J^*(cal(C)_"bi") = "SSE"(B) = w (w + 1)$. Conversely, the true 3-means optimum $cal(C)^* = { {-L (times w)}, {-1, +1}, {+L (times w)} }$ attains $J^*(cal(C)^*) = (-1)^2 + (+1)^2 = 2$. Thus
  $ J^*(cal(C)_"bi") / J^*(cal(C)^*) = (w (w + 1)) / 2 > w^2 / 2 >= M. $
  Finally, the leaf centroids of $cal(C)_"bi"$ are $mu_1 = -(w + 2)$, $mu_2 = -1$, and $mu_3 = w + 1$. For $x = +1 in C_3$, $|x - mu_2| = 2 < w = |x - mu_3|$ for all $w >= 3$, proving that $x = +1$ violates Voronoi optimality.
]

#remark(title: "Value of Global Post-Refinement")[
  In Theorem 2, running a *single* global Lloyd refinement step (Algorithm 1, Step 8) immediately reassigns $x = +1$ to $mu_2 = -1$, jumping directly to the global optimum $cal(C)^*$ and collapsing cost from $w(w + 1)$ to $2$.
]

#subsec("C", "Height-Parameterized Complexity Analysis")

#theorem(title: "Exact Tree-Height Complexity of Bisecting K-Means", origin: "Original Proof")[
  Let $H = max_(v in "leaves"(cal(T))) "depth"(v)$ be the height of the binary tree $cal(T)$ built by Algorithm 1 with $R$ trials of at most $T$ iterations per split. Then:
  1. The total arithmetic complexity is $O(n d R T H + k log k)$.
  2. For balanced trees ($H = O(log k)$), runtime is $O(n d R T log k)$, achieving a strict asymptotic speedup over flat $k$-means $O(n k d T)$ whenever $R log k = o(k)$.
  3. For maximally skewed comb trees ($H = k - 1$), runtime degrades to $Theta(n k d R T)$ for $k <= n / 2$.
]

#proof[
  Let $"int"(cal(T))$ be the $k - 1$ internal nodes of $cal(T)$. At every depth $ell in {0, 1, dots, H - 1}$, the clusters ${C_u : u in "int"(cal(T)), "depth"(u) = ell}$ are pairwise disjoint subsets of ${1, dots, n}$, so
  $ sum_(u in "int"(cal(T))) |C_u| = sum_(ell=0)^(H-1) sum_("depth"(u) = ell) |C_u| <= sum_(ell=0)^(H-1) n = n H. $
  Each internal node $u$ runs $R$ trials of 2-means (each costing $O(|C_u| dot 2 d T)$ flops) plus $O(log k)$ heap operations, summing to $O(d R T sum_u |C_u| + k log k) = O(n d R T H + k log k)$. When every split is balanced ($|A|, |B| in [c |C|, (1 - c)|C|]$ for constant $c in (0, 1/2]$), $H = O(log k)$, proving Claim 2. Conversely, if each split peels off a singleton point ($|A_m| = 1, |B_m| = n - m$), then $H = k - 1$ and $sum_(m=0)^(k-2) (n - m) = (k - 1) n - binom(k - 1, 2) = Theta(n k)$ for $k <= n / 2$, proving Claim 3.
]

#subsec("D", "Local Hardness, Hierarchical Cost, and Outlier Peeling")

#theorem(title: "Local Split Hardness and Dasgupta Cost Behavior", origin: "Classical — Cited")[
  1. *Local NP-Hardness:* Computing the globally optimal 2-means split of an internal cluster $C_u subset.eq bb(R)^d$ is NP-hard in general dimension $d$ @aloise2009np, and individual Lloyd 2-means trials can exhibit super-polynomial iteration counts in worst-case planar constructions @vattani2011kmeans.
  2. *Dasgupta Tree Cost:* Under Dasgupta's hierarchical clustering objective @dasgupta2016cost, Moseley and Wang @moseley2017approximation proved that while agglomerative average linkage achieves a constant $1/3$-approximation on the dual revenue objective, top-down Bisecting K-Means fails to admit a constant worst-case Dasgupta approximation factor across general Euclidean point clouds, despite producing more balanced clusters than standard K-Means on high-dimensional text corpora @steinbach2000comparison.
]

#proof(label: "Proof Synthesis (Cited)")[
  Aloise et al. @aloise2009np establish NP-hardness of exact 2-means via reduction from graph cut optimization in $bb(R)^(|E|)$, proving that every individual internal bisection in Algorithm 1 is already NP-hard to solve globally. Moseley and Wang @moseley2017approximation analyze hierarchical objectives by contrasting cut-balance penalties in Dasgupta's cost $sum_(i, j) w_(i j) |T(i, j)|$: because 2-means minimizes squared Euclidean spread rather than normalized similarity cut, adversarial configurations with distant sparse tails induce imbalanced top-level splits that inflate internal subtree sizes $|T(i, j)|$ on high-similarity pairs by an unbounded factor relative to the optimal tree.
]

#proposition(title: "Isometry Invariance and Outlier-Peeling Breakdown", origin: "Original Proof")[
  1. Algorithm 1 is strictly invariant under rigid Euclidean isometries $T(x) = Q x + b$ ($Q in O(d)$), yet not invariant under anisotropic diagonal scaling $S(x) = "diag"(s_1, dots, s_d) x$, which alters both 2-means hyperplanes and the max-SSE heap priority order.
  2. The centroid breakdown point is $epsilon^* = 1 / n$; however, under exact 2-means splitting with max-SSE selection, sending a single outlier $norm(x_1)_2 -> infinity$ strictly forces Step 1 to peel ${x_1}$ into an isolated singleton leaf, leaving all remaining $n - 1$ clean observations uncorrupted across the remaining $k - 1$ clusters.
]

#proof[
  (1) Because $Q in O(d)$ preserves Euclidean norms and commutes with means, every $"SSE"(C)$ and every 2-means objective is preserved identically under $T$. Under anisotropic scaling, stretching one coordinate $s_1 -> infinity$ inflates $"SSE"(C)$ selectively on clusters aligned with axis $1$, reversing both heap extraction order and hyperplane orientation.

  (2) Let $X' = cal(X) backslash {x_1}$ be fixed with finite diameter $D$ and $"SSE"(X') = S_0 < infinity$. Any 2-partition $(A, B)$ that places $x_1$ alongside at least one clean point $x_j in X'$ satisfies $"SSE"(A) + "SSE"(B) >= 1 / 2 norm(x_1 - x_j)_2^2 -> infinity$ as $norm(x_1)_2 -> infinity$. By contrast, isolating $A = {x_1}$ and $B = X'$ yields total cost $0 + S_0 = S_0 < infinity$. Hence for sufficiently large $norm(x_1)_2$, the unique optimal root split peels ${x_1}$ with diverging centroid $mu_1 = x_1$ (proving $epsilon^* = 1 / n$) while shielding $X'$ from centroid drift at the price of one cluster slot.
]

#sec("III", "Complexity, Geometry and Operational Profile")

By executing only binary subproblems on progressively smaller active subsets, Bisecting K-Means reduces per-point distance checks from $k$ to $2$ per tree level, yielding $O(n d R T log k)$ complexity on balanced hierarchies (Theorem 3) with linear $O(n d + k d)$ memory. Because each internal split runs on a localized subset $C_u$, its working set rapidly fits into L2/L3 CPU cache as depth increases, and the binary dendrogram $cal(T)$ simultaneously provides an $O(d log k)$ top-down routing index for assigning unseen test points without scanning all $k$ centroids.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Theoretical Property], [Characterization for Bisecting K-Means],
    [Time / Space Complexity], [$O(n d R T log k)$ balanced, $O(n k d R T)$ skewed; $O((n+k)d)$ space.],
    [Geometric Inductive Bias], [Recursive half-space BSP tree; leaves are convex polyhedra (not necessarily Voronoi).],
    [Hierarchy & Routing], [Produces full binary dendrogram $cal(T)$ enabling $O(d log k)$ approximate lookup.],
    [Approximation Guarantee], [Unbounded worst-case $k$-means ratio without refinement (Theorem 2); strict monotonic drop (Theorem 1).],
    [Outlier Behavior], [Breakdown $epsilon^* = 1 / n$ via singleton peeling (Proposition 1), protecting clean clusters.],
    [Recommended Upgrades], [Add global Lloyd pass (Step 8) for CVT optimality; X-Means (Alg. 10) for automatic $k$.],
  ),
  caption: [Summary of theoretical guarantees and operational profile of Bisecting K-Means.],
)

#sec("IV", "Engineering and Scientific Applications")

Bisecting K-Means thrives in high-dimensional, large-$k$, or multi-resolution applications where a full binary taxonomy is as valuable as the flat leaf partition itself.

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Large-Scale Text & Document Taxonomy Construction:* Originally validated on TREC and Reuters corpora @steinbach2000comparison, high-dimensional TF–IDF document spaces exhibit sparse orthogonal topics. Flat $k$-means frequently leaves degenerate or vastly uneven clusters, whereas max-SSE bisection systematically subdivides dominant topics first, producing balanced topic hierarchies with $O(log k)$ traversal depth @dhillon2001concept.
2. *Hierarchical Vector Quantization & Tree-Structured Codebooks:* In speech, audio, and neural embedding compression, Tree-Structured Vector Quantization (TSVQ) mirrors Bisecting K-Means directly (closely linked to LBG splitting, Alg. 30 @linde1980algorithm): storing the $k - 1$ internal split hyperplanes enables $O(d log k)$ nearest-codeword encoding instead of $O(k d)$ exhaustive search.
3. *Web-Scale Embedding Sharding & Retrieval Trees:* When $k$ reaches $10^4$–$10^5$ partitions in billion-vector search indices, running flat $k$-means is computationally prohibitive. Balanced divisive bisection builds coarse-to-fine partition trees in $O(n d log k)$ time with naturally bounded shard sizes.
4. *Multi-Resolution Geospatial & Market Regionalization:* Retail territory design and logistics zoning require nested regional hierarchies (nation $->$ regions $->$ districts). Because Lemma 1 drives each split to maximize spatial variance reduction, Bisecting K-Means produces compact, contiguous convex subregions at every granularity.
5. *Exploratory Interactive Drill-Down Analytics:* Analysts inspecting customer or gene-expression cohorts benefit from on-demand splitting of a single high-variance cluster without globally reshuffling previously interpreted sibling segments—an operational capability impossible in flat K-Means.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Bisecting K-Means degrades under four conditions: (a) *early irreversible mis-splits* on multi-modal symmetric structures (Theorem 2), requiring global Lloyd post-refinement; (b) *non-convex manifolds* (interlocking arcs or rings), because every split boundary is a linear hyperplane; (c) *heavy outlier showers*, where max-SSE selection wastes cluster budget peeling individual extreme points (Proposition 1); and (d) *unknown cluster count $k$*, where fixed $k$ must be replaced by statistical stopping criteria such as X-Means (Alg. 10) or G-Means @hamerly2003gmeans.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Why Bisecting K-Means Fits], [Caution / Alternative],
    [Document Corpora], [Balanced topic trees; $O(log k)$ depth], [Use Spherical (Alg. 12) for cosine],
    [TSVQ / Codebooks], [Fast $O(d log k)$ binary routing], [Add global Lloyd refinement pass],
    [Large-$k$ Indexing], [Sublinear in $k$ vs. flat K-Means], [Guard against skewed comb trees],
    [Spatial Zoning], [Nested convex regional drill-down], [Use DBSCAN (Alg. 31) for irregular arcs],
    [Noisy / Outlier Data], [Peels isolated outliers cleanly], [Pre-filter outliers or use Trimmed (Alg. 26)],
  ),
  caption: [Application suitability and domain engineering guide for Bisecting K-Means.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Bisecting K-Means bridges flat centroid partitioning (Category I) and agglomerative variance-minimizing hierarchies (Category III), sharing the exact same incremental inertia functional as Ward's method in reverse temporal direction.

#proposition(title: "Exact Objective Duality with Ward's Agglomerative Linkage", origin: "Original Derivation")[
  Let $Delta_"Ward" (A, B)$ denote Ward's agglomerative merge cost @ward1963hierarchical between disjoint clusters $A, B subset.eq cal(X)$. Then for every cluster $C = A union B$,
  $ Delta "SSE"(C -> A, B) = Delta_"Ward" (A, B) = (|A| |B|) / (|A| + |B|) norm(mu_A - mu_B)_2^2. $
  Consequently, greedy Bisecting K-Means (maximizing $Delta "SSE"$ top-down) and Ward's method (Alg. 59, minimizing $Delta_"Ward"$ bottom-up) are exact dual greedy heuristics over the same global WCSS functional $J^*(cal(C))$, yet are not tree-equivalent for $k >= 3$.
]

#proof[
  By Ward's definition @ward1963hierarchical, merging $A$ and $B$ increases total inertia by $"SSE"(A union B) - "SSE"(A) - "SSE"(B)$, which by Lemma 1 equals $(|A| |B| / (|A| + |B|)) norm(mu_A - mu_B)_2^2 = Delta "SSE"(C -> A, B)$ identically. Non-equivalence for $k >= 3$ follows directly from Theorem 2: on the 1D dataset $cal(X) = {-L (times w), -1, +1, +L (times w)}$ ($L = w + 2$, $w >= 3$), Ward's bottom-up linkage first merges all co-located points at $-L$ and $+L$ at zero cost and next merges ${-1}$ with ${+1}$ at cost $2 < w(w + 1)$, recovering the exact global optimum $cal(C)^*$ ($J^* = 2$), whereas top-down Bisecting K-Means is trapped in $cal(C)_"bi"$ ($J^* = w(w + 1)$).
]

This duality also links Algorithm 9 directly to LBG binary codebook splitting (Alg. 30 @linde1980algorithm), X-Means (Alg. 10 @pelleg2000xmeans, which replaces fixed $k$ with BIC-gated bisections and adds global refinement), and G-Means @hamerly2003gmeans (which gates each binary split using a projected Anderson–Darling Gaussianity test).

Empirically on the ClusterOverview 4-dataset suite, Bisecting K-Means achieves strong overall accuracy at low latency (mean runtime $0.292$ ms): Iris $"Sil" = 0.447$ and $"ARI" = 0.583$, Customers $"Sil" = 0.351$, Multi-Shape $"Sil" = 0.713$ ($"ARI" = "N/A"$), and High-D Subspace+Noise $"Sil" = 0.300$ ($"ARI" = "N/A"$). Notably, its Iris ARI ($0.583$) substantially surpasses K-Prototypes ($0.430$) and X-Means ($0.306$) because repeated multi-trial 2-means cleanly isolates the linearly separated setosa cluster at Step 1 and then stably partitions the remaining versicolor/virginica mass, while its slightly higher runtime ($0.292$ ms vs. $0.154$ ms) reflects the $R$-trial binary search overhead formalized in Theorem 3.
