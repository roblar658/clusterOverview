#import "../template/ieee.typ": *

#alg_entry(
  id: 57,
  name: "Complete Linkage Clustering (CLINK)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1977",
  paradigm: "Agglomerative Minimax-Diameter / Clique Hierarchy",
  complexity: $O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Let $(X, d)$ be a finite dissimilarity space with $X = {x_1, dots, x_n}$ and symmetric, non-negative, zero-diagonal dissimilarity $d : X times X -> bb(R)_(>= 0)$. *Complete Linkage* (also called farthest-neighbor linkage) represents the dual extreme of Single Linkage (Alg. 56) within the agglomerative hierarchy: the dissimilarity between two disjoint nonempty clusters $A, B subset.eq X$ is defined by their worst-case (maximum) pairwise separation:
$
  d_C (A, B) = max_(x in A, \, y in B) d(x, y).
$
When clusters $A$ and $B$ merge, the updated complete-linkage dissimilarity to any third active cluster $C subset.eq X without (A union B)$ obeys the Lance--Williams recurrence @lance1967general with parameters $alpha_A = alpha_B = 1 / 2$, $beta = 0$, and $gamma = +1 / 2$:
$
  d_C (A union B, C) &= 1 / 2 d_C (A, C) + 1 / 2 d_C (B, C) \
  &quad + 1 / 2 abs(d_C (A, C) - d_C (B, C)) \
  &= max { d_C (A, C), \, d_C (B, C) }.
$
Define the *diameter* of a subset $C subset.eq X$ by $"diam"(C) = max_(x, y in C) d(x, y)$ (with $"diam"({x}) = 0$). For any disjoint $A, B subset.eq X$, the merged diameter satisfies
$
  "diam"(A union B) = max { "diam"(A), \, "diam"(B), \, d_C (A, B) }.
$
Because greedy complete linkage always merges the pair $(A, B)$ minimizing $d_C (A, B)$ and merge heights never decrease (proved in Theorem 1), each step greedily minimizes the *maximum cluster diameter* among all single-merge coarsenings of the current partition.

An important algorithmic distinction separates two $O(n^2)$ constructions frequently conflated in the literature. Defays' *CLINK* algorithm @defays1977clink adapts Sibson's pointer representation $(pi, lambda)$ from SLINK to obtain an $O(n^2)$-time, $O(n)$-space streaming hierarchy; however, unlike SLINK, Defays' CLINK does *not* always produce the exact greedy complete-linkage dendrogram and depends on input order. By contrast, the *Nearest-Neighbor Chain (NN-Chain)* algorithm @mullner2011modern exploits Bruynooghe's Reducibility property to compute the *exact* complete-linkage dendrogram in optimal $Theta(n^2)$ time (using $Theta(n^2)$ dissimilarity storage).

#definition(title: "Complete-Linkage Hierarchy and Superdominant Ultrametric")[
  Given $(X, d)$, an ultrametric $u in cal(U)(X)$ is *superdominant* (or *dominating*) if $u(x, y) >= d(x, y)$ for all $x, y in X$. Denote the superdominant cone by $cal(U)_(+)(d) = { u in cal(U)(X) : u >= d }$. The *exact complete-linkage dendrogram* is the merge sequence $(A_t, B_t, h_t)_(t=1)^(n-1)$ where $(A_t, B_t) in arg min_(A != B in cal(C)_(n-t+1)) d_C (A, B)$ and $h_t = d_C (A_t, B_t)$. Its induced cophenetic dissimilarity $u_C (x, y)$ equals the merge height $h_t$ at which $x$ and $y$ first belong to the same cluster.
]

#algbox(
  title: "Exact O(n^2) Complete Linkage via Nearest-Neighbor Chain",
  inputs: [Symmetric dissimilarity matrix $D in bb(R)_(>= 0)^(n times n)$ on $X = {x_1, dots, x_n}$; tie-breaking total order $prec$; target $k$.],
  outputs: [Exact merge sequence ${(A_t, B_t, h_t)}_(t=1)^(n-1)$ and $k$-partition $cal(C)_k$.],
)[
  1. Initialize active clusters $cal(C) <- {{x_1}, dots, {x_n}}$, working matrix $d_C <- D$, empty stack $S <- []$, and $t <- 1$.
  2. *While* $abs(cal(C)) > 1$:
     - If $S$ is empty, push an arbitrary active cluster $C in cal(C)$ onto $S$.
     - Let $A$ be the top of $S$. Find $B = arg min_(C in cal(C) without {A}) d_C (A, C)$ (breaking ties deterministically via $prec$).
     - *If* $B$ is already on $S$ (necessarily immediately below $A$, forming a reciprocal nearest-neighbor pair):
       - Pop $A$ and $B$ from $S$; record merge $(A_t, B_t, h_t) <- (A, B, d_C (A, B))$.
       - Form $U <- A union B$ and update $cal(C) <- (cal(C) without {A, B}) union {U}$.
       - For every $C in cal(C) without {U}$, set $d_C (U, C) <- max { d_C (A, C), d_C (B, C) }$.
       - Increment $t <- t + 1$.
     - *Else:* push $B$ onto $S$.
  3. Sort merges by ascending $h_t$ (guaranteed non-decreasing) and stop after $n - k$ merges to output $cal(C)_k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Superdominant Ultrametric and Clique-Partition Characterization")

We first prove that exact complete linkage guarantees non-decreasing merge heights, exact diameter tracking, superdominance, and a clique-partition interpretation in threshold graphs.

#theorem(title: "Monotonicity, Superdominance, and Clique Characterization", origin: "Original Proof")[
  Let $(A_t, B_t, h_t)_(t=1)^(n-1)$ be the exact complete-linkage merge sequence on $(X, d)$, and let $u_C$ be its cophenetic dissimilarity. Then:
  1. *Reducibility and Zero Inversions:* $h_1 <= h_2 <= dots <= h_(n-1)$, and $u_C$ is a valid ultrametric on $X$.
  2. *Exact Diameter Identity:* Every cluster $U$ created at step $t$ satisfies $"diam"(U) = h_t$.
  3. *Superdominance:* $u_C in cal(U)_(+)(d)$, i.e., $u_C (x, y) >= d(x, y)$ for all $x, y in X$.
  4. *Clique Partitioning:* For any threshold $delta >= 0$, cutting the dendrogram at height $delta$ yields a partition $cal(C)(delta)$ in which every cluster $C in cal(C)(delta)$ is a clique of the threshold graph $G_delta = (X, { {x, y} : d(x, y) <= delta })$.
]

#proof[
  *Part 1:* First verify Bruynooghe's Reducibility property: for any disjoint active clusters $A, B, C$,
  $
    d_C (A union B, C) &= max { d_C (A, C), \, d_C (B, C) } \
    &>= min { d_C (A, C), \, d_C (B, C) }.
  $
  At step $t$, $(A_t, B_t)$ achieves the global minimum $h_t = min_(P != Q) d_C (P, Q)$ over active clusters. For step $t + 1$, every active pair $(P, Q)$ either already existed at step $t$ (hence $d_C (P, Q) >= h_t$) or involves the newly created cluster $U_t = A_t union B_t$ and some $C$. In the latter case, reducibility yields
  $
    d_C (U_t, C) >= min { d_C (A_t, C), \, d_C (B_t, C) } >= h_t.
  $
  Taking the minimum over all pairs at step $t + 1$ proves $h_(t+1) >= h_t$. Since merge heights are non-decreasing along every rootward path in the merge tree, for any $x, y, z in X$ the lowest common ancestor $"lca"(x, z)$ lies at or below the higher of $"lca"(x, y)$ and $"lca"(y, z)$, so $u_C (x, z) <= max { u_C (x, y), u_C (y, z) }$, confirming $u_C in cal(U)(X)$.

  *Part 2:* Induct on $t$. Base singletons have $"diam"({x}) = 0$. If $U_t = A_t union B_t$ merges clusters formed at earlier steps $s_A, s_B < t$, the inductive hypothesis and Part 1 give $"diam"(A_t) = h_(s_A) <= h_t$ and $"diam"(B_t) = h_(s_B) <= h_t$. Therefore
  $
    "diam"(U_t) &= max { "diam"(A_t), \, "diam"(B_t), \, d_C (A_t, B_t) } \
    &= max { h_(s_A), \, h_(s_B), \, h_t } = h_t.
  $

  *Parts 3 and 4:* For any distinct $x, y in X$, let $U_t$ be the cluster formed at the step where $x$ and $y$ first unite. Then $u_C (x, y) = h_t = "diam"(U_t) >= d(x, y)$, proving $u_C in cal(U)_(+)(d)$. Finally, any cluster $C in cal(C)(delta)$ was formed at some height $h <= delta$, so $d(x, y) <= "diam"(C) = h <= delta$ for all $x, y in C$; hence every pair in $C$ is connected by an edge in $G_delta$, making $C$ a clique of $G_delta$.
]

#subsec("B", "Non-Uniqueness of Minimal Superdominant Ultrametrics and Suboptimality")

In sharp contrast to Single Linkage's unique subdominant ultrametric $u_S$ (Alg. 56), the cone $cal(U)_(+)(d)$ of superdominant ultrametrics lacks a unique minimum, and greedy complete linkage can fail to achieve the minimum-diameter $k$-clustering by nearly a factor of $2$ even on four collinear points.

#theorem(title: "Non-Uniqueness in U_+(d) and Greedy Diameter Suboptimality", origin: "Original Counterexample")[
  Fix $epsilon in (0, 1 / 2)$ and consider the four-point Euclidean metric space $X = {x_1, x_2, x_3, x_4} subset bb(R)$ with coordinates $x_1 = 0$, $x_2 = 1$, $x_3 = 2 - epsilon$, and $x_4 = 3 - epsilon$. Then:
  1. There exist two incomparable minimal elements $u^((1)), u^((2)) in cal(U)_(+)(d)$ (so $cal(U)_(+)(d)$ has no unique minimum).
  2. Exact complete linkage greedily merges ${x_2, x_3}$ first at height $1 - epsilon$, forcing every resulting $2$-clustering to have maximum cluster diameter $2 - epsilon$, whereas the optimal $2$-clustering ${{x_1, x_2}, {x_3, x_4}}$ achieves maximum diameter $1$. Thus the greedy approximation ratio for minimum-diameter $2$-clustering approaches $2$ even for $n = 4$ in $bb(R)^1$.
]

#proof[
  The pairwise Euclidean distances are $d(x_1, x_2) = 1$, $d(x_2, x_3) = 1 - epsilon$, $d(x_3, x_4) = 1$, $d(x_1, x_3) = d(x_2, x_4) = 2 - epsilon$, and $d(x_1, x_4) = 3 - epsilon$.

  *Part 1:* Define $u^((1))$ by the hierarchy $({x_1, x_2} "at" 1, \, {x_3, x_4} "at" 1, \, X "at" 3 - epsilon)$ and $u^((2))$ by $({x_2, x_3} "at" 1 - epsilon, \, {x_1, x_2, x_3} "at" 2 - epsilon, \, X "at" 3 - epsilon)$. Both satisfy the strong triangle inequality and dominate $d$. Observe that $u^((1))(x_2, x_3) = 3 - epsilon > 1 - epsilon = u^((2))(x_2, x_3)$, whereas $u^((1))(x_1, x_2) = 1 < 2 - epsilon = u^((2))(x_1, x_2)$. Hence neither dominates the other. Moreover, in any $u in cal(U)_(+)(d)$ with $u <= u^((1))$, we must have $u(x_1, x_2) >= 1 = u^((1))(x_1, x_2)$, $u(x_3, x_4) >= 1 = u^((1))(x_3, x_4)$, and since $u(x_1, x_4) >= 3 - epsilon$ while $u(x_1, x_2) = 1 < 3 - epsilon$, the ultrametric inequality forces $u(x_2, x_4) = u(x_1, x_3) = u(x_2, x_3) = 3 - epsilon$; thus $u^((1))$ is minimal. An analogous check confirms $u^((2))$ is minimal.

  *Part 2:* At Step 1, the unique minimum pairwise distance is $d(x_2, x_3) = 1 - epsilon < 1$, so exact complete linkage *must* merge $A = {x_2, x_3}$. At Step 2, the active complete-linkage dissimilarities are $d_C ({x_1}, A) = max(1, 2 - epsilon) = 2 - epsilon$, $d_C (A, {x_4}) = max(2 - epsilon, 1) = 2 - epsilon$, and $d_C ({x_1}, {x_4}) = 3 - epsilon$. Regardless of tie-breaking, Step 2 merges either $x_1$ or $x_4$ into $A$, yielding a $2$-clustering of maximum diameter $2 - epsilon$. In contrast, ${{x_1, x_2}, {x_3, x_4}}$ has maximum diameter $max(1, 1) = 1$, giving ratio $(2 - epsilon) / 1 -> 2$ as $epsilon -> 0^(+)$. Note that minimizing maximum cluster diameter for general $k >= 3$ is NP-hard, and greedy complete linkage can exhibit $Omega(log k)$ approximation factors on adversarial metrics @dasgupta2016cost.
]

#subsec("C", "Exactness and Quadratic Complexity of NN-Chain vs. CLINK")

#theorem(title: "Exactness and Amortized O(n^2) Time of NN-Chain Complete Linkage", origin: "Original Proof")[
  Under deterministic tie-breaking, the NN-Chain algorithm (Algorithm Box) constructs the exact complete-linkage dendrogram in $Theta(n^2)$ time and $Theta(n^2)$ space @mullner2011modern. Conversely, Defays' $O(n^2)$-time, $O(n)$-space CLINK pointer algorithm @defays1977clink is not permutation-invariant and does not in general match exact complete linkage.
]

#proof[
  *Exactness:* Suppose $(A, B)$ are reciprocal nearest neighbors (RNNs) under $d_C$, meaning $d_C (A, B) <= d_C (A, C)$ and $d_C (A, B) <= d_C (B, C)$ for all $C in cal(C) without {A, B}$ (with strictness enforced by the tie-breaking order $prec$). We claim that if any other disjoint pair $(P, Q)$ merges before $(A, B)$ in global-minimum order, $(A, B)$ remains an RNN pair with identical merge height $d_C (A, B)$. Indeed, after merging $P$ and $Q$ into $U = P union Q$, reducibility (Theorem 1) gives
  $
    d_C (A, U) = max { d_C (A, P), \, d_C (A, Q) } >= d_C (A, B),
  $
  and identically $d_C (B, U) >= d_C (A, B)$. Thus merging $(P, Q)$ never introduces a cluster strictly closer to $A$ or $B$ than each other, nor alters $d_C (A, B)$. By induction, the RNN pair $(A, B)$ is guaranteed to merge at exact height $d_C (A, B)$ in the global dendrogram; merging RNNs immediately upon discovery only reorders independent subtrees at their true heights.

  *Termination and $Theta(n^2)$ Time:* Along the stack $S = (C_1, dots, C_r)$, edge dissimilarities are non-increasing and strictly decreasing under the deterministic tie-breaker $prec$. Because $cal(C)$ is finite and $prec$ breaks ties asymmetrically, no cycle of length $>= 3$ can occur, so the chain necessarily terminates at an RNN pair $(C_(r-1), C_r)$. Moreover, popping $C_(r-1)$ and $C_r$ upon merging them into $U = C_(r-1) union C_r$ preserves the validity of the remaining stack prefix $(C_1, dots, C_(r-2))$: by reducibility, for every $j <= r - 3$,
  $
    d_C (C_j, U) &>= min { d_C (C_j, C_(r-1)), \, d_C (C_j, C_r) } \
    &>= d_C (C_j, C_(j+1)),
  $
  so $C_(j+1)$ remains the exact nearest neighbor of $C_j$. Let $P$ denote the total number of stack pushes across the entire execution. Each of the $n - 1$ merges pops exactly $2$ clusters from $S$, and at termination $abs(S_"final") = 1$ (or $0$ if emptied prior to loop exit), giving
  $
    P - 2(n - 1) <= 1 ==> P <= 2n - 1.
  $
  Since each push scans at most $abs(cal(C)) <= n$ active distances and each of the $n - 1$ merges performs $O(abs(cal(C)))$ Lance--Williams updates, the total number of distance operations is bounded by $(3n - 2) n = Theta(n^2)$.

  *CLINK Non-Exactness:* In Defays' CLINK @defays1977clink, each incoming point $x_(m+1)$ updates the pointer representation greedily using only the existing $m$-point clusters. On the four-point line $X = {0, 1, 2 - epsilon, 3 - epsilon}$ from Theorem 2 presented in order $(x_1, x_2, x_4, x_3)$, after three points CLINK merges ${x_1, x_2}$ at height $1$; upon arrival of $x_3 = 2 - epsilon$, $x_3$ has distance $1$ to $x_4$ and $2 - epsilon$ to ${x_1, x_2}$, causing CLINK to form ${{x_1, x_2}, {x_3, x_4}}$ at height $1$ --- strictly disagreeing with exact complete linkage's merge of ${x_2, x_3}$ at $1 - epsilon$.
]

#subsec("D", "Monotone Equivariance and Extreme-Outlier Perturbation")

#proposition(title: "Ordinal Equivariance and Non-Lipschitz Outlier Amplification", origin: "Original Proof")[
  1. *Ordinal equivariance:* For any strictly increasing bijection $phi : [0, infinity) -> [0, infinity)$ with $phi(0) = 0$ on a tie-free dissimilarity space $(X, d)$, exact complete linkage satisfies $u_C (phi compose d) = phi compose u_C (d)$.
  2. *Outlier amplification:* Suppose cluster $A$ of diameter $delta$ absorbs a single outlier $z$ at distance $R > 2 delta$ from $A$. Then for every other cluster $B$ within distance $D < R - delta$ of $A$, the complete-linkage distance jumps by at least $d_C (A union {z}, B) - d_C (A, B) >= R - 2 D - delta$. Moreover, unlike single linkage, $u_C$ is *discontinuous* (not Lipschitz) in the Gromov--Hausdorff metric @carlsson2010characterization.
]

#proof[
  *Part 1:* Because $phi$ is strictly increasing, $max_(x in A, y in B) phi(d(x, y)) = phi(max_(x in A, y in B) d(x, y))$. By induction on merge step $t$, all pairwise comparisons between active complete-linkage dissimilarities under $phi compose d$ preserve exact order relative to $d$, so identical pairs merge at heights $phi(h_t)$, giving $u_C (phi compose d) = phi compose u_C (d)$.

  *Part 2:* Choose $x in A$ and $y in B$ realizing $d(x, y) <= D$. By the triangle inequality (in a metric space),
  $
    d_C (A union {z}, B) >= d(z, y) >= d(z, x) - d(x, y) >= R - D,
  $
  whereas $d_C (A, B) <= d(x, y) + "diam"(A) + "diam"(B)$. When $B$ is a singleton ($"diam"(B) = 0$), $d_C (A, B) <= D + delta$, so $d_C (A union {z}, B) - d_C (A, B) >= R - 2 D - delta$, which diverges linearly with $R$. For Gromov--Hausdorff discontinuity, perturb the four-point line ${0, 1, 2, 3}$ by $+- eta$ to ${0, 1, 2 - eta, 3 - eta}$ versus ${0, 1 - eta, 2, 3 - eta}$: the Gromov--Hausdorff distance between the two inputs is $O(eta) -> 0$, yet the complete-linkage ultrametrics jump by $1 - O(eta)$ as $eta -> 0^(+)$ because the initial tie between ${x_1, x_2}$ and ${x_2, x_3}$ breaks oppositely @carlsson2010characterization.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational and Memory Profile:* Exact complete linkage via NN-Chain runs in $Theta(n^2)$ time and $Theta(n^2)$ space (requiring the active pairwise dissimilarity matrix to execute $max$ updates). Generic heap-based implementations incur $O(n^2 log n)$ time and higher constant factors due to stale heap entries and tie-resolution bookkeeping. When memory is constrained ($n > 50{,}000$), Defays' CLINK @defays1977clink reduces space to $Theta(n)$ at the cost of order-dependent heuristic approximation (Theorem 3).

*Geometric Inductive Bias:* Complete linkage enforces *maximum-diameter compactness*: every cluster at threshold $delta$ is guaranteed to be a clique of $G_delta$. Consequently, it completely suppresses single-linkage chaining and produces balanced, compact, spherical-like groupings in Euclidean space. However, its extremal maximum operator makes it intolerant of non-convex geometries (severing crescents, rings, and elongated filaments into diameter-bounded segments) and vulnerable to distant outliers that inflate a cluster's effective radius once absorbed.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$Theta(n^2)$ exact via NN-Chain @mullner2011modern; $O(n^2 log n)$ via binary heaps.],
    [Space Complexity], [$Theta(n^2)$ for exact complete linkage; $Theta(n)$ for Defays' approximate CLINK.],
    [Determinism], [Deterministic given tie-breaking; sensitive to exact ties (non-Lipschitz).],
    [Invariances], [Exact invariance to strictly increasing transforms and positive scaling (tie-free).],
    [Cluster Geometry], [Compact, bounded-diameter, roughly spherical cliques; breaks elongated shapes.],
    [Monotonicity], [Strictly inversion-free; yields a superdominant ultrametric $u_C >= d$.],
    [Outlier Sensitivity], [Resistant to bridging chains; severe radius distortion once an outlier is absorbed.],
    [Primary Bottleneck], [$Theta(n^2)$ dissimilarity matrix storage and tie/heap maintenance overhead.],
  ),
  caption: [Operational and theoretical profile of Complete Linkage Clustering (CLINK).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Operational Taxonomic Unit (OTU) Clustering in Microbiome Sequencing:* Standard 16S rRNA amplicon pipelines define species-level OTUs via a strict $97\%$ pairwise sequence identity threshold ($d <= 0.03$). Because Theorem 1 guarantees every cluster at cut $delta = 0.03$ is a clique of $G_(0.03)$, complete linkage enforces the biological requirement that *all* sequences within an OTU share at least $97\%$ identity, avoiding the chimeric mega-OTUs produced by single linkage.
2. *Cheminformatics and Pharmacophore Diversity Binning:* In high-throughput virtual screening, compound libraries are partitioned under Tanimoto dissimilarity $d_T = 1 - T_c$ so that a single representative molecule can stand in for its entire chemical bin. Complete linkage's strict diameter bound $"diam"(C) <= delta$ certifies that every unselected compound lies within guaranteed Tanimoto radius $delta$ of its bin prototype.
3. *Wireless Sensor Network Bounded-Diameter Cell Partitioning:* Multi-hop wireless ad-hoc networks require cluster heads with strict worst-case intra-cluster propagation delay or transmission range. Complete linkage guarantees that every pair of nodes in a cluster communicates within the maximum range threshold $delta$ without transitive hop blow-up.
4. *Near-Duplicate Document and Plagiarism Grouping:* When clustering news articles or legal filings under Jaccard or cosine dissimilarity, transitive chaining via boiler-plate text merges wholly unrelated documents. Complete linkage guarantees uniform all-pairs lexical similarity inside each deduplication bucket.
5. *Worst-Case ($L_infinity$) Vector Quantization and Color Palette Compression:* When maximum perceptual color error (rather than average squared error) governs image quality, complete linkage in CIELAB space directly bounds the worst-case color distortion $"diam"(C)$ within every quantized palette cell.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Complete linkage underperforms on curved or multi-scale manifolds (where a true cluster's geodesic length exceeds the inter-cluster gap, forcing premature fragmentation), on datasets containing wide dynamic ranges of cluster diameters (large diffuse clusters are broken into small pieces matching dense clusters' diameters), and in the presence of unpruned extreme outliers.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [16S rRNA OTU Cliques], [Exact all-pairs identity bound via Theorem 1], [Order/tie sensitivity; use UPGMA (Alg. 58)],
    [Chemical Diversity Bins], [Hard Tanimoto diameter guarantee per bin], [Prototype optimality; use Minimax @bien2011minimax],
    [Document Deduplication], [Eliminates transitive boilerplate chaining], [Memory for $n > 10^5$; use LSH pre-blocking],
    [Curved / Ring Manifolds], [Poor fit; clique bias fractures non-convex shapes], [Prefer SLINK (Alg. 56) or HDBSCAN\*],
    [Multi-Scale Gaussian Blobs], [Large-variance clusters get prematurely split], [Prefer Ward (Alg. 59) or GMM @fraley2002mclust],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Complete Linkage.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Complete linkage anchors the upper boundary of the Lance--Williams ultrametric envelope, standing in formal duality to Single Linkage (Alg. 56) and motivating prototype-centered generalizations such as Minimax Linkage @bien2011minimax.

#proposition(title: "Ultrametric Sandwich Duality and Minimax Linkage Interpolation", origin: "Original Derivation")[
  For any finite dissimilarity space $(X, d)$, let $u_S$, $u_U$, and $u_C$ denote the cophenetic ultrametrics of Single (Alg. 56), Average (Alg. 58), and Complete (Alg. 57) linkage. Then:
  1. *Ultrametric sandwich:* $u_S (x, y) <= d(x, y) <= u_C (x, y)$ for all $x, y in X$, with equality $u_S = d = u_C$ if and only if $d$ is already an ultrametric (in which case $u_S = u_U = u_C = d$).
  2. *Minimax radius bracket:* Define the minimax radius $r(C) = min_(p in C) max_(x in C) d(x, p)$ used by Minimax Linkage @bien2011minimax. On any metric space $(X, d)$, every cluster $C subset.eq X$ satisfies
     $
       1 / 2 "diam"(C) <= r(C) <= "diam"(C),
     $
     showing that Minimax Linkage maintains complete linkage's bounded-diameter guarantee within a factor of $2$ while identifying an explicit cluster prototype $p in C$.
]

#proof[
  *Part 1:* Theorem 1 of Alg. 56 established $u_S <= d$, and Theorem 1 above established $d <= u_C$, immediately yielding $u_S <= d <= u_C$. If $d$ is an ultrametric, then $d in cal(U)_(-)(d)$ forces $u_S >= d$ (by maximality of $u_S$), so $u_S = d$; similarly, at every step of complete linkage on an ultrametric, all cross-distances between two distinct ultrametric balls are identical, so $d_C (A, B) = d_U (A, B) = d_S (A, B)$ and $u_C = u_U = u_S = d$. Conversely, if $u_S = u_C$, then $u_S = d = u_C in cal(U)(X)$, so $d$ is an ultrametric.

  *Part 2:* For any prototype $p in C$, $max_(x in C) d(x, p) <= max_(x, y in C) d(x, y) = "diam"(C)$, so taking the minimum over $p in C$ gives $r(C) <= "diam"(C)$. Conversely, let $p^(star) in C$ attain $r(C)$ and let $x, y in C$ attain $"diam"(C) = d(x, y)$. By the triangle inequality,
  $
    "diam"(C) = d(x, y) <= d(x, p^(star)) + d(p^(star), y) <= 2 r(C),
  $
  proving $1 / 2 "diam"(C) <= r(C) <= "diam"(C)$.
]

Complete linkage also links naturally to correlation clustering, graph clique-cover partitioning, and the remaining agglomerative benchmarks: Average Linkage (Alg. 58, which relaxes $max$ to arithmetic mean) and Ward's Method (Alg. 59, which replaces diameter minimization by Euclidean variance minimization).

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized four-dataset suite, Complete Linkage records *Iris* $"Sil" = 0.444$, $"ARI" = 0.570$; *Customers* $"Sil" = 0.323$; *Multi-Shape* $"Sil" = 0.420$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.257$ ($"ARI" = "N/A"$); and a mean runtime of *5.391 ms*. These numbers vividly exhibit the dual trade-off proved in Section II. On *Customers*, complete linkage's anti-chaining clique bias rescues Silhouette to $0.323$ — nearly a $5 times$ gain over SLINK's collapsed $0.067$ — closely tracking UPGMA ($0.334$) and Ward ($0.347$). On *Iris*, suppressing inter-species bridges nudges $"ARI"$ slightly above SLINK and UPGMA ($0.570$ vs. $0.560$), though enforcing compact diameter boundaries reduces Silhouette ($0.444$ vs. $0.489$) and still trails Ward's variance-optimal $"ARI"$ of $0.696$. Conversely, on *Multi-Shape*, complete linkage records the lowest Silhouette of all four hierarchical methods ($0.420$, well below UPGMA's $0.452$, SLINK's $0.609$, and Ward's $0.712$), directly reflecting the fracture of elongated non-convex geometries into bounded-diameter cliques (Theorem 1). On *High-D*, it achieves $0.257$, virtually tied with SLINK ($0.254$) and slightly ahead of UPGMA ($0.245$), while trailing Ward ($0.292$). Finally, its mean runtime of $5.391$ ms is more than double that of SLINK ($2.343$ ms), UPGMA ($2.549$ ms), and Ward ($2.344$ ms), reflecting higher dissimilarity-matrix update and heap/tie-maintenance overhead in exact complete linkage.
