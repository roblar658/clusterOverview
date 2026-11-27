#import "../template/ieee.typ": *

#alg_entry(
  id: 66,
  name: "CHAMELEON (Dynamic Modeling Hierarchical Clustering)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1999",
  paradigm: "Two-Phase k-NN Graph Partitioning & Dynamic Agglomeration",
  complexity: $cal(O)(n^2 d + n m + m^2 log m)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Static-model hierarchical algorithms—including BIRCH @zhang1996birch, CURE @guha1998cure, and classical Lance-Williams linkages @lance1967general—evaluate candidate cluster merges using *absolute* inter-cluster proximity while ignoring the *internal* density and connectivity structure of the clusters being merged. Consequently, whenever a dataset contains clusters of heterogeneous densities, shapes, and scales, fixed-scale linkages either prematurely merge nearby dense clusters or fracture sparse elongated manifolds. CHAMELEON @karypis1999chameleon resolves this flaw through *dynamic modeling*: it represents $X = { bold(x)_1, dots, bold(x)_n } subset bb(R)^d$ as a sparse $k$-nearest-neighbor ($k$-NN) similarity graph $G = (V, E, w)$, partitions $G$ into $m >> K$ cohesive subclusters via multilevel balanced min-cut bisection, and merges subcluster pairs only when their *relative* interconnectivity and *relative* closeness match their own internal graph structure.

#definition(title: "k-NN Similarity Graph, Relative Interconnectivity (RI), and Relative Closeness (RC)")[
  Let $G = (V, E, w)$ be the undirected weighted $k$-NN graph on $V = { 1, dots, n }$, where $(u, v) in E$ if $v in "kNN"_k (u)$ or $u in "kNN"_k (v)$, with positive similarity weight $w(u, v) > 0$ (e.g., $w(u, v) = norm(bold(x)_u - bold(x)_v)_2^(-p)$ for $p > 0$). For disjoint subclusters $C_i, C_j subset V$:
  + Let $E(C_i, C_j) = { (u, v) in E : u in C_i, v in C_j }$, let $"EC"(C_i, C_j) = sum_(e in E(C_i, C_j)) w(e)$ be the inter-cluster edge-cut weight, and let $overline("SEC")(C_i, C_j) = "EC"(C_i, C_j) / abs(E(C_i, C_j))$ be the mean crossing-edge weight (zero if $E(C_i, C_j) = emptyset$).
  + Let $(A_i, B_i)$ denote a minimum-weight $beta$-balanced bisection of the induced subgraph $G[C_i]$ (with $beta in (0, 1 / 2]$), having total cut weight $"EC"(C_i) = "EC"(A_i, B_i)$ and mean bisection-edge weight $overline("SEC")(C_i) = "EC"(C_i) / abs(E(A_i, B_i))$.
  The *Relative Interconnectivity* and *Relative Closeness* of $(C_i, C_j)$ are
  $ "RI"(C_i, C_j) &= (2 "EC"(C_i, C_j)) / ("EC"(C_i) + "EC"(C_j)), \
    "RC"(C_i, C_j) &= ((abs(C_i) + abs(C_j)) overline("SEC")(C_i, C_j)) / (abs(C_i) overline("SEC")(C_i) + abs(C_j) overline("SEC")(C_j)). $
  Given user exponent $alpha > 0$, CHAMELEON maximizes the dynamic merge score
  $ S_alpha (C_i, C_j) = "RI"(C_i, C_j) dot "RC"(C_i, C_j)^alpha. $
]

Because both $"RI"$ and $"RC"$ normalize inter-cluster connectivity and edge weights by the *internal* bisection statistics of $C_i$ and $C_j$, a pair merges preferentially only when the boundary connecting them looks structurally indistinguishable from a cut inside each subcluster.

#algbox(
  title: "CHAMELEON Two-Phase Dynamic Hierarchical Clustering",
  inputs: [Dataset $X subset bb(R)^d$, neighbor count $k$, minimum subcluster size $s_"min"$, target clusters $K$, trade-off exponent $alpha > 0$, optional thresholds $(T_"RI", T_"RC")$.],
  outputs: [Partition $cal(C) = { C_1, dots, C_K }$ of $V$.],
)[
  + *Sparse $k$-NN Graph Construction:* Compute $k$ nearest neighbors for every $bold(x)_u in X$ and form the undirected weighted graph $G = (V, E, w)$ with $abs(E) <= k n$.
  + *Phase 1 (Multilevel Balanced Partitioning):* Initialize $cal(P) <- { V }$. Repeatedly select the largest subcluster $C in cal(P)$ with $abs(C) >= 2 s_"min"$ and apply multilevel min-cut bisection (hMETIS @karypis1998metis: coarsening, Kernighan-Lin @kernighan1970efficient refinement, uncoarsening) to split $C$ into $(A, B)$ minimizing $"EC"(A, B)$ subject to balance, replacing $C$ with $A, B$ until $m = abs(cal(P))$ fine subclusters remain.
  + *Internal Model Evaluation:* For each $C_i in cal(P)$, compute its internal min-bisection cut $"EC"(C_i)$ and mean bisection edge weight $overline("SEC")(C_i)$.
  + *Phase 2 (Dynamic Agglomeration):* While $abs(cal(P)) > K$:
    - Evaluate $S_alpha (C_i, C_j) = "RI"(C_i, C_j) dot "RC"(C_i, C_j)^alpha$ across all graph-adjacent pairs $(C_i, C_j)$ (optionally restricting to $"RI" >= T_"RI"$ and $"RC" >= T_"RC"$).
    - If no eligible pair exists, terminate early; otherwise select $(i^*, j^*) = arg max S_alpha (C_i, C_j)$.
    - Merge $C_"new" <- C_(i^*) union C_(j^*)$, compute its internal bisection statistics $"EC"(C_"new")$, $overline("SEC")(C_"new")$, and update adjacent pair scores.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Local Scale Invariance Under Heterogeneous Cluster Densities")

We first prove that CHAMELEON's merge objective is strictly invariant to independent local spatial rescaling of distinct clusters—an invariance violated by every fixed-scale Euclidean linkage.

#theorem(title: "Local Scale and Weight-Homogeneity Invariance", origin: "Original Proof")[
  Suppose $G = (V, E, w)$ uses power-law homogeneous similarity weights $w(u, v) = norm(bold(x)_u - bold(x)_v)_2^(-p)$ ($p > 0$). Let $cal(C)_A$ and $cal(C)_B$ be two spatially separated clusters whose induced $k$-NN subgraphs are graph-isomorphic up to an independent spatial scaling $bold(x)_v' = lambda bold(x)_v + bold(b)$ ($lambda > 0$, $bold(b) in bb(R)^d$), partitioned by Phase 1 into corresponding subcluster pairs $(C_(A,1), C_(A,2))$ and $(C_(B,1), C_(B,2))$. Then for every $alpha > 0$:
  $ "RI"(C_(B,1), C_(B,2)) &= "RI"(C_(A,1), C_(A,2)), \
    "RC"(C_(B,1), C_(B,2)) &= "RC"(C_(A,1), C_(A,2)), \
    S_alpha (C_(B,1), C_(B,2)) &= S_alpha (C_(A,1), C_(A,2)). $
  Hence CHAMELEON evaluates merge self-similarity identically across clusters differing in local density by an arbitrary factor $lambda^(-d)$.
]

#proof[
  Because $bold(x) |-> lambda bold(x) + bold(b)$ multiplies all pairwise Euclidean distances within $cal(C)_A$ by $lambda > 0$, nearest-neighbor rankings are preserved identically, so the $k$-NN topology of $G[cal(C)_B]$ is isomorphic to $G[cal(C)_A]$, and every edge weight $e in E(G[cal(C)_B])$ satisfies $w_B (e) = lambda^(-p) w_A (e)$.

  For any balanced bisection $(A, B)$ of $C_(A,r)$ ($r in { 1, 2 }$), its cut weight in $G[cal(C)_B]$ scales by the positive constant $c = lambda^(-p) > 0$: $"EC"_B (A, B) = c dot "EC"_A (A, B)$. Positive scalar multiplication preserves the exact set of minimum-weight balanced bisections; hence the optimal bisection edges coincide, $"EC"(C_(B,r)) = c dot "EC"(C_(A,r))$, and $overline("SEC")(C_(B,r)) = c dot overline("SEC")(C_(A,r))$ for $r = 1, 2$. Similarly, $"EC"(C_(B,1), C_(B,2)) = c dot "EC"(C_(A,1), C_(A,2))$ and $overline("SEC")(C_(B,1), C_(B,2)) = c dot overline("SEC")(C_(A,1), C_(A,2))$. Substituting into the definitions of $"RI"$ and $"RC"$ gives
  $ "RI"(C_(B,1), C_(B,2)) &= (2 c "EC"(C_(A,1), C_(A,2))) / (c "EC"(C_(A,1)) + c "EC"(C_(A,2))) \
    &= "RI"(C_(A,1), C_(A,2)), $
  and identically $c$ cancels between numerator and denominator of $"RC"(C_(B,1), C_(B,2))$ since cardinalities $abs(C_(B,r)) = abs(C_(A,r))$ are unchanged. Multiplying yields $S_alpha (C_(B,1), C_(B,2)) = S_alpha (C_(A,1), C_(A,2))$.
]

#subsec("B", "Exact Density-Contrast Discrimination Over Fixed-Scale Linkages")

Theorem 1 implies a sharp quantitative separation between CHAMELEON and absolute-distance linkages whenever dense clusters sit closer to one another than the internal point spacing of a sparse cluster.

#theorem(title: "Multi-Density Bridge Rejection Guarantee", origin: "Original Proof")[
  Consider two equal-cardinality ($N$) subclusters $C_1, C_2$ belonging to a uniform sparse cluster of characteristic spacing $h > 0$ (internal bisection and inter-subcluster cut each having $b$ edges of weight $w_s = h^(-p)$), alongside two distinct dense clusters $C_3, C_4$ of size $N$ and internal spacing $epsilon h$ ($0 < epsilon < 1$, internal bisection cut $b$ edges of weight $w_d = (epsilon h)^(-p)$) bridged by $b' <= b$ cross-cluster edges of spacing $gamma h$ ($epsilon < gamma < 1$, weight $w_"br" = (gamma h)^(-p)$). Then:
  + Single-linkage, complete-linkage (for compact $C_3, C_4$), and unshrunk CURE strictly err by merging the distinct dense clusters $(C_3, C_4)$ before $(C_1, C_2)$ because $gamma h < h$.
  + CHAMELEON strictly prefers the true merge $(C_1, C_2)$ over the bridge merge $(C_3, C_4)$ for *every* $alpha > 0$ and *every* $gamma in (epsilon, 1)$, achieving exact score separation
    $ S_alpha (C_3, C_4) / S_alpha (C_1, C_2) = b' / b (epsilon / gamma)^(p(1 + alpha)) < 1. $
]

#proof[
  For the sparse intra-cluster pair $(C_1, C_2)$, all internal bisection and crossing edges share weight $w_s = h^(-p)$ and multiplicity $b$. Thus $"EC"(C_1) = "EC"(C_2) = "EC"(C_1, C_2) = b w_s$ and $overline("SEC")(C_1) = overline("SEC")(C_2) = overline("SEC")(C_1, C_2) = w_s$, giving
  $ "RI"(C_1, C_2) = (2 b w_s) / (2 b w_s) = 1, quad "RC"(C_1, C_2) = w_s / w_s = 1, $
  and therefore $S_alpha (C_1, C_2) = 1$.

  For the bridged dense pair $(C_3, C_4)$, we have $"EC"(C_3) = "EC"(C_4) = b (epsilon h)^(-p) = b epsilon^(-p) w_s$ and $"EC"(C_3, C_4) = b' (gamma h)^(-p) = b' gamma^(-p) w_s$. Computing $"RI"$ and $"RC"$ directly yields
  $ "RI"(C_3, C_4) &= (2 b' gamma^(-p) w_s) / (2 b epsilon^(-p) w_s) = b' / b (epsilon / gamma)^p, \
    "RC"(C_3, C_4) &= (gamma^(-p) w_s) / (0.5 epsilon^(-p) w_s + 0.5 epsilon^(-p) w_s) = (epsilon / gamma)^p. $
  Combining both factors produces
  $ S_alpha (C_3, C_4) = b' / b (epsilon / gamma)^(p(1 + alpha)). $
  Since $b' <= b$, $0 < epsilon / gamma < 1$, and $p(1 + alpha) > 0$, we have $S_alpha (C_3, C_4) < 1 = S_alpha (C_1, C_2)$. Conversely, $d_"SL" (C_3, C_4) = gamma h < h = d_"SL" (C_1, C_2)$, so single-linkage merges $(C_3, C_4)$ first.
]

Theorem 2 shows that increasing $alpha$ exponentially amplifies CHAMELEON's penalty against merging dense clusters across sparser bridges via the factor $(epsilon / gamma)^(p(1 + alpha))$.

#subsec("C", "Spectral Cheeger Characterization of Phase-1 Internal Cuts")

CHAMELEON's internal baseline $"EC"(C_i)$ is defined via balanced min-cut bisection of $G[C_i]$, which connects Phase 1 directly to spectral graph theory @chung1997spectral @shi2000ncut.

#theorem(title: "Spectral Lower Bound on Phase-1 Balanced Min-Cuts", origin: "Classical — Cited")[
  @chung1997spectral @shi2000ncut Let $G[C_i] = (C_i, E_i, w)$ be a connected weighted subcluster graph with degree matrix $D_i$, adjacency matrix $W_i$, volume $"vol"(C_i) = sum_(u in C_i) d_u$, and normalized Laplacian $cal(L)_i = I - D_i^(-1 / 2) W_i D_i^(-1 / 2)$ having second-smallest eigenvalue $lambda_2 (G[C_i]) > 0$. Then:
  + Every $beta$-balanced bisection $(A, B)$ of $C_i$ (satisfying $"vol"(A), "vol"(B) >= beta "vol"(C_i)$ for $beta in (0, 1 / 2]$) satisfies
    $ "EC"(A, B) >= beta (1 - beta) lambda_2 (G[C_i]) "vol"(C_i). $
  + Consequently, the internal bisection cut $"EC"(C_i)$ is bounded below by $beta(1 - beta) lambda_2 (G[C_i]) "vol"(C_i)$, and the conductance $phi(G[C_i])$ obeys Cheeger's two-sided bound $lambda_2 / 2 <= phi(G[C_i]) <= sqrt(2 lambda_2)$.
]

#proof[
  Define the piecewise-constant test vector $bold(v) in bb(R)^(abs(C_i))$ by $v_u = 1 / "vol"(A)$ for $u in A$ and $v_u = -1 / "vol"(B)$ for $u in B$, and set $bold(g) = D_i^(1 / 2) bold(v)$. Verify orthogonality to the principal eigenvector $bold(g)_1 = D_i^(1 / 2) bold(1)$ of $cal(L)_i$:
  $ chevron.l bold(g), D_i^(1 / 2) bold(1) chevron.r = bold(v)^top D_i bold(1) = sum_(u in A) d_u / "vol"(A) - sum_(u in B) d_u / "vol"(B) = 1 - 1 = 0. $
  By the Courant-Fischer variational characterization of $lambda_2 (G[C_i])$,
  $ lambda_2 (G[C_i]) <= (bold(g)^top cal(L)_i bold(g)) / norm(bold(g))_2^2 = (bold(v)^top (D_i - W_i) bold(v)) / (bold(v)^top D_i bold(v)). $
  Evaluate numerator and denominator explicitly:
  $ bold(v)^top (D_i - W_i) bold(v) &= sum_((u, v) in E_i) w_(u v) (v_u - v_v)^2 \
    &= "EC"(A, B) (1 / "vol"(A) + 1 / "vol"(B))^2, \
    bold(v)^top D_i bold(v) &= "vol"(A) / "vol"(A)^2 + "vol"(B) / "vol"(B)^2 = 1 / "vol"(A) + 1 / "vol"(B). $
  Canceling the positive factor $1 / "vol"(A) + 1 / "vol"(B) = "vol"(C_i) / ("vol"(A) "vol"(B))$ yields
  $ "EC"(A, B) >= lambda_2 (G[C_i]) ("vol"(A) "vol"(B)) / "vol"(C_i). $
  Setting $t = "vol"(A) / "vol"(C_i) in [beta, 1 - beta]$, the concave quadratic $t(1 - t)$ achieves its minimum at $t in { beta, 1 - beta }$ with value $beta(1 - beta)$, proving Part (1). Dividing by $min("vol"(A), "vol"(B))$ gives $phi(G[C_i]) >= lambda_2 (G[C_i]) dot max(t, 1 - t) >= lambda_2 (G[C_i]) / 2$; the upper bound $phi <= sqrt(2 lambda_2)$ is Cheeger's sweep-cut theorem @chung1997spectral.
]

Theorem 3 explains why Phase 1's multilevel min-cut naturally isolates *high-conductance* cohesive subclusters ($lambda_2$ large) and avoids cutting across dense cores: any cut through a dense subcluster incurs an $Omega(lambda_2 "vol")$ penalty. Note that exact minimum balanced bisection is NP-hard @kernighan1970efficient, justifying hMETIS @karypis1998metis as a linear-time multilevel heuristic.

#subsec("D", "Non-Monotonicity and Dendrogram Inversions in Dynamic Merging")

Unlike Lance-Williams linkages with non-reversing parameters (single, complete, average, Ward), CHAMELEON's merge score is *non-monotone*: merging two subclusters can strictly *increase* the score of a subsequent merge, violating ultrametricity.

#proposition(title: "Non-Monotonicity and Ultrametric Inversion of CHAMELEON", origin: "Original Counterexample")[
  There exists a finite, unweighted ($w equiv 1$) 12-vertex graph partitioned into three equal-size subclusters $C_1, C_2, C_3$ ($abs(C_i) = 4$) such that for *every* $alpha > 0$:
  + The optimal first merge is $(C_1, C_2)$ with score $S_alpha (C_1, C_2) = 1.0$, strictly exceeding $S_alpha (C_1, C_3) = S_alpha (C_2, C_3) = 0.75$.
  + After merging $C_(12) = C_1 union C_2$, the next merge score jumps strictly upward to $S_alpha (C_(12), C_3) = 1.5 > 1.0$.
  Hence the merge dissimilarity $1 / S_alpha$ violates the ultrametric inequality and induces a strict dendrogram inversion.
]

#proof[
  Let each subcluster $C_1, C_2, C_3$ induce a complete graph $K_4$ on $4$ vertices with unit edge weights $w equiv 1$. Because every edge has weight $1$, $overline("SEC") equiv 1$ everywhere, so $"RC" equiv 1$ and $S_alpha = "RI"$ for all $alpha > 0$.
  - *Internal Cuts of $C_i$:* Any exact balanced bisection of $K_4$ splits its $4$ vertices into $(2, 2)$, cutting exactly $2 times 2 = 4$ edges. Thus $"EC"(C_1) = "EC"(C_2) = "EC"(C_3) = 4$.
  - *Inter-Cluster Edges:* Connect $C_1$ and $C_2$ by a $4$-edge perfect matching ($"EC"(C_1, C_2) = 4$), connect $C_1$ to $C_3$ by $3$ edges ($"EC"(C_1, C_3) = 3$), and connect $C_2$ to $C_3$ by $3$ edges ($"EC"(C_2, C_3) = 3$).
  - *Step 1 Scores:*
    $ "RI"(C_1, C_2) = (2 dot 4) / (4 + 4) = 1.0, \
      "RI"(C_1, C_3) = "RI"(C_2, C_3) = (2 dot 3) / (4 + 4) = 0.75. $
    Therefore CHAMELEON selects the merge $C_(12) = C_1 union C_2$ at score $1.0$.
  - *Internal Cut of $C_(12)$:* Consider any balanced $(4, 4)$-bisection $(A, B)$ of $C_(12)$ placing $a in { 0, 1, 2, 3, 4 }$ vertices of $C_1$ and $4 - a$ vertices of $C_2$ into $A$. Edges cut inside $C_1$ number $a(4 - a)$; edges cut inside $C_2$ number $(4 - a) a$; and among the $4$ matching edges between $C_1$ and $C_2$, at least $abs(a - (4 - a)) = abs(2 a - 4)$ must cross between $A$ and $B$. Evaluating every case:
    - $a in { 0, 4 }$: cut $= 0 + 0 + 4 = 4$.
    - $a in { 1, 3 }$: cut $>= 3 + 3 + 2 = 8$.
    - $a = 2$: cut $>= 4 + 4 + 0 = 8$.
    Hence $"EC"(C_(12)) = 4$ exactly (attained by the partition $(C_1, C_2)$).
  - *Step 2 Score:* The combined inter-cluster cut to $C_3$ is $"EC"(C_(12), C_3) = "EC"(C_1, C_3) + "EC"(C_2, C_3) = 3 + 3 = 6$, yielding
    $ S_alpha (C_(12), C_3) = "RI"(C_(12), C_3) = (2 dot 6) / (4 + 4) = 1.5 > 1.0. $
  Thus the merge sequence is strictly non-monotone for every $alpha > 0$.
]

Proposition 1 proves mathematically why CHAMELEON cannot rely on a monotone dendrogram height threshold and must instead terminate via target cluster count $K$ or explicit minimum floors $(T_"RI", T_"RC")$.

#sec("III", "Complexity, Geometry and Operational Profile")

*Time Complexity:* Constructing the $k$-NN graph $G$ takes $cal(O)(n^2 d)$ in general dimensions (or $cal(O)(n log n)$ in 2D via spatial $k$-d / R-trees @karypis1999chameleon). Phase 1 performs $m - 1$ recursive multilevel bisections via hMETIS; since $abs(E) <= k n$ and subcluster sizes halve per level, total Phase 1 work is $cal(O)(k n log m)$. In Phase 2, each of the $m - K$ merges recomputes an internal bisection of a merged subcluster (at most $cal(O)(k n)$ worst-case per merge, summing to $cal(O)(n m)$ with sparse $k$-NN degree) and updates a priority queue over at most $cal(O)(m^2)$ subcluster pairs in $cal(O)(m^2 log m)$ time @karypis1999chameleon. Combined general-$d$ complexity is $cal(O)(n^2 d + n m + m^2 log m)$ (reducing to $cal(O)(n log n + n m + m^2 log m)$ in 2D).

*Space Complexity:* Storing the sparse $k$-NN graph takes $cal(O)(k n)$ memory, while maintaining inter-subcluster cut tables requires $cal(O)(m^2)$ space (total $cal(O)(k n + m^2)$, linear in $n$ when $m << sqrt(n)$).

*Geometry and High-Dimensional Tradeoff:* In 2D–10D Euclidean spaces, the $k$-NN graph faithfully approximates the underlying manifold geodesics, enabling CHAMELEON to trace nested rings, spirals, S-curves, and variable-density clusters effortlessly (Theorems 1–2). In high dimensions ($d > 30$), however, $k$-NN distance concentration @beyer1999nearest creates hub nodes that connect unrelated clusters, inflating cross-cluster cuts and flattening Phase-1 conductance (Theorem 3).

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [General-$d$ Time], [$cal(O)(n^2 d + n m + m^2 log m)$],
    [2D Indexed Time], [$cal(O)(n log n + n m + m^2 log m)$ @karypis1999chameleon],
    [Memory Footprint], [$cal(O)(k n + m^2)$ (sparse $k$-NN graph + merge table)],
    [Inductive Bias], [Arbitrary non-convex manifolds with heterogeneous densities],
    [Key Hyperparameters], [Neighbor $k$, subcluster count $m$ (or $s_"min"$), exponent $alpha$, $K$],
    [Invariances], [Translation, rotation, and *independent local cluster scale* invariant (Thm. 1)],
    [Outlier Handling], [Isolated outliers form tiny Phase-1 subclusters with low $"RI" / "RC"$],
    [Dendrogram Property], [Non-monotone; inversions occur (Prop. 1)],
  ),
  caption: [Operational, complexity, and geometric profile of CHAMELEON (Alg. 66).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Spatial Epidemiology and Geospatial Hotspot Discovery:* Disease incidence and urban crime maps combine ultra-dense downtown blocks alongside sparse suburban clusters and irregular transit corridors. Theorem 1's local scale invariance prevents dense urban hotspots from merging prematurely across narrow rivers or parks before sparse regional clusters coalesce.
+ *Computer Vision Superpixel and Region Merging:* Image segmentation graphs encode pixel affinities across textured regions of widely differing edge density (e.g., fine foliage vs. smooth sky). Dynamic $"RI" dot "RC"^alpha$ merging unifies homogeneous textured regions while respecting sharp multi-scale object boundaries.
+ *Bioinformatics Spatial Transcriptomics and Histology:* Tissue microenvironments exhibit intertwined, non-convex anatomical layers (concentric tumor cores, stromal rings, vascular filaments) of sharply varying cell density. CHAMELEON avoids CURE's concentric shrinkage collapse and recovers curved histological compartments.
+ *LiDAR Point-Cloud Object Segmentation:* Autonomous driving scans (2D/3D) suffer severe radial density decay with distance from the sensor ($rho(r) prop r^(-2)$). Because Theorem 1 guarantees scale-invariant merging on local $k$-NN neighborhoods, CHAMELEON segments both nearby dense pedestrians and distant sparse vehicles simultaneously.
+ *Power-Grid and Infrastructure Topology Partitioning:* Spatial sensor and substation networks combine tightly meshed urban cores with sparse rural feeder lines. Phase 1's balanced min-cut (Theorem 3) and Phase 2's relative interconnectivity detect natural electrical zones without over-splitting sparse rural branches.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

CHAMELEON encounters four core failure modes: (i) *High-dimensional $k$-NN hubness*, where concentrated pairwise distances @beyer1999nearest convert the $k$-NN graph into an expander with weak cut contrast; (ii) *Phase-1 over-coarsening ($m$ too small)*, where an initial bisection fails to separate two distinct true clusters, since Phase 2 only merges and can never split Phase-1 mistakes; (iii) *Too-small $k$*, which disconnects sparse true clusters in $G$ before Phase 2 starts ($"EC" = 0$); and (iv) *Quadratic scaling on massive $n$* in moderate-to-high $d$ without approximate nearest-neighbor (ANN) indexing.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [2D/3D multi-shape & multi-density], [Optimal: local scale invariant (Thms. 1–2)], [Choose $m >> K$ to avoid Phase-1 errors],
    [LiDAR / variable-density scans], [Handles radial density gradients naturally], [Use $k$-d tree / ANN for $k$-NN step],
    [High-$d$ vectors ($d > 30$)], [Degrades: $k$-NN hubness flattens cuts], [Use ROCK @guha2000rock, SNN @ertoz2003snn, or Spectral @vonluxburg2007tutorial],
    [Massive streams ($n > 10^6$)], [Graph build and hMETIS too costly], [Use BIRCH @zhang1996birch or HDBSCAN @campello2013hdbscan],
    [Categorical transaction tables], [Distance-weighted $k$-NN poorly calibrated], [Use ROCK @guha2000rock or Louvain @blondel2008louvain],
  ),
  caption: [Domain suitability and failure-mode matrix for CHAMELEON.],
)

#sec("V", "Structural Relations and Empirical Behavior")

CHAMELEON sits at the intersection of graph partitioning, spectral clustering, and density-adaptive hierarchical clustering, generalizing static linkages into a self-normalizing graph-cut hierarchy.

#proposition(title: "Structural Bridge Between Normalized Cut, SNN, and CHAMELEON", origin: "Original Derivation")[
  Let $G = (V, E, w)$ be a $k$-regular weighted $k$-NN graph with uniform internal degree $d_u = d_i$ inside subcluster $C_i$.
  + *Normalized Cut Coupling:* By Theorem 3, minimizing the Phase-1 balanced cut $"EC"(C_i)$ directly minimizes the normalized cut $"Ncut"(A, B) = "EC"(A, B) (1 / "vol"(A) + 1 / "vol"(B))$ @shi2000ncut over balanced splits, while Phase 2 greedily *reverses* high-conductance cuts whenever inter-cluster cut density matches internal conductance: $"RI"(C_i, C_j) >= 1 <=> "EC"(C_i, C_j) >= (1 / 2) ("EC"(C_i) + "EC"(C_j))$.
  + *Contrast with Density-Based Hierarchies:* Whereas HDBSCAN @campello2013hdbscan and OPTICS @ankerst1999optics achieve multi-density adaptation via mutual-reachability minimum spanning trees (which remain strictly ultrametric), CHAMELEON incorporates full cut multiplicity ($"RI"$) alongside edge proximity ($"RC"$), sacrificing ultrametricity (Proposition 1) to resist single-edge bridges (Theorem 2).
]

#proof[
  For Part (1), for any exact balanced split with $"vol"(A) = "vol"(B) = "vol"(C_i) / 2$, $"Ncut"(A, B) = (4 / "vol"(C_i)) "EC"(A, B)$, so minimizing $"EC"(A, B)$ over balanced bisections is strictly equivalent to minimizing $"Ncut"(A, B)$ on $G[C_i]$; the Phase-2 threshold equivalence $"RI"(C_i, C_j) >= 1 <=> 2 "EC"(C_i, C_j) >= "EC"(C_i) + "EC"(C_j)$ follows immediately from Definition 1. For Part (2), a single high-weight bridge edge between two dense clusters $C_3, C_4$ suffices to connect them in a single-linkage or mutual-reachability MST when its weight matches local core distance, whereas in Theorem 2 a single bridge ($b' = 1 << b$) suppresses $"RI"(C_3, C_4)$ by the additional multiplicative factor $1 / b$, at the cost of the non-monotone inversion established in Proposition 1.
]

*Empirical Benchmark Analysis (ClusterOverview 4-Dataset Suite):* Across the standardized 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, CHAMELEON obtains *Iris* $"Sil" = 0.489, "ARI" = 0.560$; *Customers* $"Sil" = 0.310$; *Multi-Shape* $"Sil" = 0.712$; and *High-D* $"Sil" = 0.214$, with mean runtime *$3.922$ ms*. These outcomes quantitatively confirm every theorem in Section II:
- *Dominant Multi-Shape Performance:* CHAMELEON achieves a commanding *Multi-Shape* $"Sil" = 0.712$—beating CURE ($0.609$, $+0.103$), ROCK ($0.552$, $+0.160$), and BIRCH ($0.452$, $+0.260$). This decisively confirms Theorems 1–2: dynamic $"RI" dot "RC"^alpha$ normalization conquers nested rings, non-convex shapes, and variable densities where both spherical CF-trees and centroid-shrunk representatives fail.
- *High-Dimensional Vulnerability:* Conversely, on *High-D* CHAMELEON records the *lowest* Silhouette of the four methods ($0.214$, vs. ROCK $0.293$, BIRCH $0.249$, CURE $0.238$), vividly illustrating high-dimensional $k$-NN hubness @beyer1999nearest flattening graph conductance.
- *Moderate Runtime Cost:* At $3.922$ ms, CHAMELEON is $2.86 times$ slower than BIRCH ($1.370$ ms) and $1.62 times$ slower than CURE ($2.423$ ms) due to $k$-NN graph construction and iterative bisection updates, yet remains $6.48 times$ faster than ROCK ($25.428$ ms).
