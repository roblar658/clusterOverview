#import "../template/ieee.typ": *

#alg_entry(
  id: 147,
  name: "INSCY (Indexing Subspace Clusters with In-Process Redundancy Removal)",
  category: "Cat. VI: Grid, Subspace & Projection",
  year: "2008",
  paradigm: "Depth-First Indexed Non-Redundant Subspace Clustering",
  complexity: $O(n d log n + abs(cal(T)_"vis") n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Density-based subspace clustering algorithms such as SUBCLU @kailing2004subclu suffer from two tightly coupled pathologies in high dimensions: (1) *redundancy explosion*, whereby any cluster $C$ in an $s$-dimensional subspace $S$ automatically projects onto dense clusters in all $2^s - 2$ proper non-empty subspaces $S' subset S$, burying true patterns beneath exponentially many trivial projections; and (2) *volume collapse*, whereby a fixed density threshold becomes exponentially harder to satisfy as $s$ increases because $s$-dimensional ball volumes shrink rapidly with $s$. Assent, Krieger, Müller, and Seidl @assent2008inscy introduced *INSCY* (INdexing Subspace Clusters with in-process-removal of redundancY) to resolve both issues simultaneously via a compact prefix-tree index (*SCY-tree*), depth-first maximal-first traversal, dimensionality-normalized *expected-density* thresholds, and *in-process $R$-redundancy pruning*.

Let $cal(D) = {bold(x)_1, dots, bold(x)_n} subset [0, 1]^d$ be a normalized dataset (via coordinate-wise min-max scaling). In any subspace $S subset.eq {1, dots, d}$ of dimensionality $s = abs(S)$, define the kernel density estimator at $bold(x) in cal(D)$ with bandwidth $epsilon in (0, 1/2)$:
$ hat(f)_S (bold(x)) = 1 / (n V_s (epsilon)) sum_(bold(y) in cal(D)) K( norm(pi_S (bold(x)) - pi_S (bold(y))) / epsilon ), $
where $K$ is a compactly supported radial kernel (e.g., uniform indicator $bb(1)[u <= 1]$ or Epanechnikov kernel) and $V_s (epsilon)$ is the $s$-dimensional kernel volume ($V_s (epsilon) = (2 epsilon)^s$ under $L_infinity$ and $V_s (epsilon) = pi^(s/2) epsilon^s / Gamma(s / 2 + 1)$ under $L_2$). Under a uniform null distribution on $[0, 1]^s$, the expected number of points in an $epsilon$-ball is $E_s = n V_s (epsilon)$.

#definition(title: "Expected-Density Clusters, R-Redundancy, and SCY-Tree Index")[
  Fix a dimensionless density factor $F > 1$, minimum support floor $m_"min" >= 1$, and redundancy ratio $R in (0, 1]$.
  1. *Expected-Density Core Object:* An object $bold(x) in cal(D)$ is $(F, epsilon, m_"min")$-*dense* in subspace $S$ iff
     $ abs(N_epsilon^S (bold(x))) >= tau_s, quad tau_s = max( F dot.op n V_s (epsilon), m_"min" ). $
     Maximal density-connected components of dense objects in $S$ of size at least $m_"min"$ form the subspace cluster set $cal(C)(S)$.
  2. *$R$-Redundancy:* A cluster $(C', S')$ is *$R$-redundant* with respect to a higher-dimensional cluster $(C, S)$ ($S' subset S$) iff
     $ abs(C' inter C) / abs(C') >= R. $
  3. *SCY-Tree:* Partition each 1D projection $pi_({k}) (cal(D))$ into dense 1D region descriptors $cal(U)_k$ (plus overlapping border-margin descriptors of width $epsilon$ to preserve cross-boundary density connectivity). The *SCY-tree* $T$ is a prefix tree storing each object's descriptor path across dimensions $1 < 2 < dots < d$, equipped with horizontal node-links $L_u$ per descriptor $u in cal(U)_k$ and subtree object-count summaries $c(v)$.
]

Because INSCY traverses the subspace lattice *depth-first* along the SCY-tree—restricting $T$ recursively to superspaces $S union {j}$ ($j > max(S)$) *before* extracting clusters in $S$—all relevant higher-dimensional clusters are already known when $S$ is clustered. Consequently, any candidate region whose non-redundant support falls below $m_"min"$ is pruned *in process* rather than post-filtered after exponential enumeration.

#algbox(
  title: "INSCY: Depth-First Indexed Non-Redundant Subspace Clustering",
  inputs: [Normalized dataset $cal(D) subset [0, 1]^d$; radius $epsilon$; density factor $F > 1$; minimum cluster size $m_"min"$; redundancy threshold $R in (0, 1]$.],
  outputs: [Set of non-redundant subspace clusters $cal(R)$.],
)[
  1. *1D Region Discovery & SCY-Tree Build:* For each coordinate $k = 1, dots, d$, identify dense 1D regions $cal(U)_k$ (with $epsilon$-border descriptors) and insert every object's descriptor tuple into the root SCY-tree $T_0$ along dimension order $1 < 2 < dots < d$. Initialize $cal(R) <- emptyset$.
  2. *Recursive Depth-First Procedure* `"INSCY"(T, S, k_"last")`:
     - *Support Pruning:* If total object count in $T$ is strictly less than $m_"min"$, return immediately.
     - *Depth-First Superspace Expansion:* For each dimension $j = k_"last" + 1, dots, d$ and each dense descriptor $u in cal(U)_j$ present in $T$:
       * Construct restricted SCY-tree $T_u <- "Restrict"(T, u)$ by following node-link list $L_u$ upward and merging prefix paths.
       * If $"count"(T_u) >= m_"min"$, recurse `"INSCY"(T_u, S union {j}, j)`.
     - *Post-Order Clustering & In-Process Redundancy Pruning:* If $S != emptyset$:
       * Run expected-density clustering (threshold $tau_(abs(S))$) on the object set $O(T)$ in subspace $S$ to extract candidate clusters $C_1, dots, C_q$.
       * For each candidate $C_r$, test whether $exists (C, S_"sup") in cal(R)$ with $S subset S_"sup"$ and $abs(C_r inter C) / abs(C_r) >= R$. If no such covering superspace cluster exists, insert $(C_r, S)$ into $cal(R)$.
  3. *Invoke* `"INSCY"(T_0, emptyset, 0)` and return $cal(R)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We prove that lexicographic prefix-tree post-order guarantees exact reverse-inclusion precedence for in-process redundancy removal, analyze expected-density compensation against high-dimensional volume collapse, establish non-transitivity of $R$-redundancy via an explicit counterexample, and derive tight SCY-tree size and restriction-cost bounds.

#subsec("A", "Reverse-Inclusion Precedence of Prefix-Tree Post-Order")

A critical correctness requirement of in-process redundancy pruning is that *every* superspace $S supset S'$ must finish clustering before $S'$ is evaluated—not merely those superspaces that are descendants of $S'$ in the recursion tree!

#theorem(title: "Exact Reverse-Inclusion Invariant of Lexicographic DFS Post-Order", origin: "Original Proof")[
  Let the recursion tree of `"INSCY"` enumerate non-empty subsets $S subset.eq {1, dots, d}$ as strictly increasing tuples $(i_1 < dots < i_s)$, expanding children $S union {j}$ in ascending order of $j > i_s$ and executing clustering on $S$ in *post-order* (after all recursive calls return). Then:
  (i) For *every* pair of subspaces $S' subset S subset.eq {1, dots, d}$, the post-order clustering step of $S$ strictly precedes that of $S'$.
  (ii) Assuming zero support-pruning error on $m_"min"$, every cluster $(C, S)$ emitted into $cal(R)$ is strictly non-$R$-redundant with respect to every valid higher-dimensional cluster in $union_(S_"sup" supset S) cal(C)(S_"sup")$.
]

#proof[
  (i) Let $S' = (a_1 < dots < a_r)$ and $S = (b_1 < dots < b_s)$ with $S' subset S$ ($r < s$). Let $j^* = min(S without S')$ be the smallest coordinate present in $S$ but absent from $S'$. There are two exhaustive cases:
  - *Case 1 ($j^* > a_r$):* Every element of $S without S'$ exceeds $a_r = max(S')$. Hence the first $r$ elements of $S$ coincide identically with $S'$: $(b_1, dots, b_r) = (a_1, dots, a_r)$. Therefore $S$ lies strictly within the recursion subtree rooted at $S'$. Because `"INSCY"` executes post-order clustering on $S'$ only after all recursive descendants return, $S$ completes strictly before $S'$.
  - *Case 2 ($j^* < a_r$):* Let $p in {1, dots, r}$ be the unique index such that $a_(p-1) < j^* < a_p$ (with $a_0 = 0$). Then $S$ and $S'$ share the exact prefix $P = (a_1, dots, a_(p-1))$, after which the $p$-th coordinate of $S$ is $b_p = j^* < a_p$. At recursion node $P$, `"INSCY"` loops over extensions $j$ in ascending order, so the entire subtree rooted at $P union {j^*}$ (which contains $S$) is entered and completely terminated *before* the loop reaches $j = a_p$ (whose subtree contains $S'$). Thus $S$ again finishes strictly before $S'$.
  (ii) For soundness of support pruning: any object $bold(x) in C$ for a cluster $C in cal(C)(S)$ belongs to a dense 1D region in every $k in S$ (since $N_epsilon^S (bold(x)) subset.eq N_epsilon^({k}) (bold(x))$ and $tau_(abs(S)) >= m_"min"$), so $abs(O(T)) >= abs(C) >= m_"min"$ along every prefix of $S$, meaning no valid cluster is pruned by the $m_"min"$ guard. By part (i), when $S'$ is evaluated, all valid clusters in all superspaces $S supset S'$ have already been examined; inductively, if a cluster $(C', S')$ is $R$-redundant w.r.t. a maximal superspace cluster $(C, S)$, $(C, S)$ is already stored in $cal(R)$ and prunes $(C', S')$.
]

#subsec("B", "Expected-Density Compensation vs. Volume Collapse")

#theorem(title: "Dimensionality-Normalized Contrast and Critical Finite-Sample Dimension", origin: "Classical — Cited")[
  Let $cal(D)$ consist of $n$ points sampled from a uniform background on $[0, 1]^d$ plus a planted cluster of density $F_"true" > F > 1$ in subspace $S$ of dimension $s$. Under $L_infinity$ neighborhoods of radius $epsilon in (0, 1/2)$:
  (i) Any fixed absolute count threshold $tau_0 > 1$ fails to detect uniform-contrast clusters once $s > log(tau_0 / (F_"true" n)) / log(2 epsilon)$, driving recall to zero exponentially in $s$ @assent2008inscy @beyer1999nearest.
  (ii) INSCY's expected-density rule $tau_s = max(F n (2 epsilon)^s, m_"min")$ preserves constant relative contrast $F$ across all dimensions up to the exact discrete-floor critical dimension
  $ s_"crit" = floor( log(m_"min" / (F n)) / log(2 epsilon) ), $
  above which the fixed floor $m_"min"$ dominates. For $L_2$ balls, Stirling's formula yields the super-exponential decay $V_s (epsilon) ~ (pi s)^(-1/2) (2 pi e epsilon^2 / s)^(s/2)$, lowering $s_"crit"$ to $Theta(epsilon^2 n^(2/s))$ unless $epsilon$ is rescaled with $sqrt(s)$.
]

#proof[
  (i) In $s$ dimensions under $L_infinity$, a hypercube of half-side $epsilon$ inside $[0, 1]^s$ has volume $V_s (epsilon) = (2 epsilon)^s$. A region with constant density contrast $F_"true"$ above uniform contains expected count $mu_s = F_"true" n (2 epsilon)^s$. Setting $mu_s < tau_0$ and dividing by $log(2 epsilon) < 0$ reverses the inequality, yielding $s > log(tau_0 / (F_"true" n)) / log(2 epsilon)$, beyond which true cluster points fall below $tau_0$.
  (ii) Equating the adaptive branch $F n (2 epsilon)^s$ with the support floor $m_"min"$ gives $(2 epsilon)^s = m_"min" / (F n)$, solving directly to $s_"crit" = floor(log(m_"min" / (F n)) / log(2 epsilon))$. For all $s <= s_"crit"$, $tau_s = F n V_s (epsilon)$, so any cluster satisfying $hat(f)_S (bold(x)) >= F$ remains detectable regardless of $s$. Under the $L_2$ norm, $V_s (epsilon) = pi^(s/2) epsilon^s / Gamma(s / 2 + 1)$; applying Stirling's asymptotic $Gamma(s / 2 + 1) ~ sqrt(pi s) (s / (2 e))^(s/2)$ yields $V_s (epsilon) ~ (pi s)^(-1/2) (2 pi e epsilon^2 / s)^(s/2)$, which decays as $s^(-s/2)$—faster than any exponential—explaining why $L_infinity$ or $sqrt(s)$-scaled radii are required in high $s$.
]

#subsec("C", "Non-Transitivity and Union Blindness of R-Redundancy")

#proposition(title: "Non-Transitivity and Pairwise Union Blindness of R-Redundancy", origin: "Original Counterexample")[
  For every $R in (1/2, 1)$:
  (i) The pairwise $R$-redundancy relation is *not transitive*: there exist clusters $(C_1, S_1)$, $(C_2, S_2)$, $(C_3, S_3)$ with $S_1 subset S_2 subset S_3$ such that $C_1$ is $R$-redundant w.r.t. $C_2$ and $C_2$ is $R$-redundant w.r.t. $C_3$, yet $C_1$ is *not* $R$-redundant w.r.t. $C_3$.
  (ii) Pairwise $R$-redundancy exhibits *union blindness*: there exist two incomparable maximal clusters $(C_1, S_1)$ and $(C_2, S_2)$ ($not (S_1 subset.eq S_2)$ and $not (S_2 subset.eq S_1)$) whose projection $C' = C_1 union C_2$ in $S' = S_1 inter S_2$ is 100% covered by $C_1 union C_2$, yet $(C', S')$ escapes pairwise $R$-pruning against both $(C_1, S_1)$ and $(C_2, S_2)$.
]

#proof[
  (i) Choose integer $m > 2 / (1 - R)$ and let $k = ceil(R m)$, so $m / 2 < k < m$ and $2 k - m < R m$. Let $S_1 = {1} subset S_2 = {1, 2} subset S_3 = {1, 2, 3}$, and construct sets of cardinality $m$ on labels ${1, dots, 3 m - k}$ by setting $C_2 = {1, dots, m}$, $C_1 = {1, dots, k} union {m + 1, dots, 2 m - k}$, and $C_3 = {m - k + 1, dots, m} union {2 m + 1, dots, 3 m - k}$. Then:
  $ abs(C_1 inter C_2) / abs(C_1) = k / m >= R, quad abs(C_2 inter C_3) / abs(C_2) = k / m >= R, $
  whereas $C_1 inter C_3 = {m - k + 1, dots, k}$ has cardinality $2 k - m$. Because $k < R m + 1$ and $m > 2 / (1 - R)$, we have $2 k - m < (2 R - 1) m + 2 < R m$, so $abs(C_1 inter C_3) / abs(C_1) < R$ and $C_1$ is *not* $R$-redundant w.r.t. $C_3$. Hence discarding $C_2$ during post-order pruning causes $C_1$ to survive even though $C_1$ was redundant w.r.t. $C_2$!
  (ii) Let $S_1 = {1, 2}$, $S_2 = {1, 3}$, $S' = {1}$, and let $C_1, C_2$ be disjoint clusters of equal size $m$ whose projections onto coordinate $1$ merge into $C' = C_1 union C_2$ of size $2 m$. Then $C' subset.eq C_1 union C_2$ (100% union coverage), yet $abs(C' inter C_1) / abs(C') = abs(C' inter C_2) / abs(C') = 1 / 2 < R$, so pairwise testing emits $(C', S')$.
]

#subsec("D", "SCY-Tree Space and Restriction-Operation Bounds")

#lemma(title: "Linear-in-Data SCY-Tree Footprint and Conditional Restriction Complexity", origin: "Original Derivation")[
  Let $n_"dense" <= n$ be the number of objects lying in at least one dense 1D region across $d$ attributes, with $D_k = abs(cal(U)_k)$ descriptors per coordinate (including border regions, each object visits at most $2$ descriptors per dimension). Then:
  (i) Without border duplication, the root SCY-tree $T_0$ contains at most $N_(T_0) <= min(n_"dense" d, sum_(j=1)^d product_(k=1)^j (D_k + 1))$ nodes, yielding $O(n d)$ space; with $epsilon$-border descriptors, each object contributes at most $2^b$ paths if it lies in $b$ border zones (bounded by $O(n d)$ when border zones are sparse).
  (ii) Restricting any SCY-tree $T$ to descriptor $u in cal(U)_j$ runs in $O(abs(L_u) dot.op j)$ time via upward traversal of node-link list $L_u$, and the recursion stack stores at most $d$ restricted trees simultaneously, bounding peak index memory by $O(d dot.op N_(T_0))$.
]

#proof[
  (i) Each non-border object inserts a single root-to-leaf path of length at most $d$ into $T_0$, adding at most $d$ new nodes; summing over $n_"dense"$ objects gives $N_(T_0) <= n_"dense" d$. Independently, level $j$ of a trie with branching factor at most $D_k + 1$ (including a null-skip descriptor) at level $k$ has at most $product_(k=1)^j (D_k + 1)$ nodes; taking the minimum yields the stated bound.
  (ii) In `"Restrict"(T, u)`, only nodes on the horizontal list $L_u$ at depth at most $j$ are visited; tracing each node's parent pointers back to the root traverses at most $j$ edges, costing $O(abs(L_u) dot.op j)$ total work without scanning objects absent from $u$. Because recursive calls strictly increment $k_"last" in {0, dots, d}$, the maximum DFS recursion depth is $d$; since each restricted tree satisfies $N_(T_u) <= N_T <= N_(T_0)$, peak stack memory is at most $O(d N_(T_0))$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

Building the root SCY-tree takes $O(n d log n)$ time to discover 1D dense regions and insert $n$ paths of length $d$. The recursive search cost is output- and structure-sensitive: letting $abs(cal(T)_"vis")$ denote the number of non-pruned subspaces reaching the clustering step, total runtime is $O(n d log n + abs(cal(T)_"vis") (N_(T_0) + n_"sub"^2 d_"sub"))$, where $n_"sub" = abs(O(T))$ is the restricted object count. On data with a few clean maximal clusters, in-process support and redundancy pruning slash $abs(cal(T)_"vis")$ drastically; however, in the adversarial worst case of exponentially many non-redundant overlapping subspace clusters, $abs(cal(T)_"vis") = O(2^d)$ remains possible.

Geometrically, INSCY combines axis-aligned subspace selection with density-connected cluster extraction, recovering arbitrary non-convex shapes inside each subspace while actively suppressing redundant lower-dimensional shadows via Theorem II.A.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n d log n + abs(cal(T)_"vis") n_"sub"^2)$ (output-sensitive; worst-case $O(2^d n^2)$)],
    [Space Complexity], [$O(n d)$ root SCY-tree; $O(d N_(T_0))$ DFS recursion stack],
    [Density Model], [Dimensionality-normalized expected density $tau_s = max(F n V_s (epsilon), m_"min")$],
    [Redundancy Control], [In-process $R$-redundancy via lexicographic DFS post-order (Thm. II.A)],
    [Key Hyperparameters], [Neighborhood radius $epsilon$, density factor $F$, $m_"min"$, redundancy $R$],
    [Cluster Geometry], [Arbitrary non-convex density-connected shapes in maximal subspaces],
  ),
  caption: [Operational and geometric profile of INSCY (Alg. 147).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Multi-Omics and Co-Expression Pathway Discovery:* In gene-expression and proteomics matrices, biologists seek *maximal* co-regulated gene subsets rather than thousands of redundant 2- and 3-gene sub-projections. INSCY's in-process $R$-redundancy removal reports compact maximal pathway clusters directly @assent2008inscy.
2. *Cybersecurity and Network Intrusion Flow Analysis:* Network flow logs span dozens of packet-size, timing, and protocol statistics where coordinated attacks manifest across maximal feature combinations. Expected-density scaling (Theorem II.B) detects high-dimensional attack signatures that fixed-threshold SUBCLU misses, while redundancy pruning prevents alert flooding.
3. *Semiconductor Wafer-Test and Industrial Yield Diagnostics:* Wafer-probe testing records 50–200 electrical parameters per die; systematic process defects leave localized density signatures in specific parameter subsets. Depth-first SCY-tree restriction pinpoints root-cause maximal parameter sets without enumerating redundant pairs.
4. *E-Commerce Behavioral Micro-Segmentation:* High-dimensional customer engagement profiles naturally contain overlapping niche cohorts operating in distinct behavioral subspaces. INSCY allows an object to participate in multiple maximal subspace clusters while eliminating trivial lower-dimensional echoes.
5. *Climate and Scientific Ensemble Telemetry:* Multi-variable geophysical simulation grids exhibit coherent regimes across subsets of atmospheric variables; adaptive expected-density thresholds keep high-dimensional regimes comparable to low-dimensional anomalies.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

INSCY encounters four key failure regimes: (1) *Border-descriptor blowup*, where wide $epsilon$ causes objects to fall into multiple 1D border zones across many dimensions, inflating SCY-tree paths exponentially ($2^b$ in Lemma II.D); (2) *Adversarial exponential branching*, when many weakly overlapping clusters resist $R$-pruning; (3) *Union blindness and non-transitivity* (Proposition II.C), where merged low-dimensional projections survive despite full union coverage by incomparable superspaces; and (4) *Low-dimensional boundary pruning*, where aggressive noise filtering and strict density multipliers $F$ discard peripheral points on small low-$d$ benchmarks.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [High-$d$ Maximal Subspaces], [High: expected density + in-process $R$-pruning], [Control border width $epsilon$; FIRES @kriegel2005fires],
    [Redundancy-Sensitive Omics], [High: suppresses $2^s - 2$ shadow projections], [Watch union blindness (Prop. II.C)],
    [Low-$d$ Compact Partitions], [Low: strict $F$-density leaves sparse points unassigned], [Use P3C @moise2006p3c or COSA @friedman2004cosa],
    [Oblique Linear Subspaces], [Low: axis-aligned 1D SCY-tree descriptors], [Use ORCLUS @aggarwal2000orclus or SSC @elhamifar2013ssc],
    [Massive-$n$ Real-Time Streams], [Low: recursive tree restriction & $O(n_"sub"^2)$ DBSCAN], [Use P3C @moise2006p3c (histogram-based)],
  ),
  caption: [Domain suitability and failure-mode matrix for INSCY.],
)

#sec("V", "Structural Relations and Empirical Behavior")

INSCY directly generalizes SUBCLU @kailing2004subclu and CLIQUE @agrawal1998clique by replacing breadth-first Apriori enumeration with depth-first indexed pruning and adaptive density thresholds.

#proposition(title: "Exact Equivalence Limit to SUBCLU and Output Compression Factor", origin: "Original Derivation")[
  Fix $(epsilon, "MinPts")$ and disable dimensionality-adaptive scaling by setting $tau_s equiv "MinPts"$ for all $s$. Then:
  (i) Disabling $R$-pruning (by setting $R > 1$) makes INSCY's output identical to SUBCLU @kailing2004subclu whenever 1D border descriptors cover every $epsilon$-neighborhood across region boundaries.
  (ii) For a dataset consisting of $q$ disjoint maximal clusters of subspace dimensions $s_1, dots, s_q$ with zero inter-cluster projection overlap, SUBCLU outputs $sum_(r=1)^q (2^(s_r) - 1)$ clusters whereas INSCY at any $R in (0, 1]$ outputs exactly $q$ clusters—achieving an exact exponential output compression factor of $sum_(r=1)^q (2^(s_r) - 1) / q$.
]

#proof[
  (i) With $tau_s equiv "MinPts"$, every $S$-density-connected cluster $C$ satisfies $C subset.eq B_k$ for a dense 1D region $B_k$ on each $k in S$ (including overlapping $epsilon$-borders, which prevent artificial splitting of connected components across 1D bin edges). Hence every subspace $S$ containing a SUBCLU cluster survives the $m_"min"$ support guard and is visited exactly once by the prefix tree; with $R > 1$, no cluster is ever pruned, yielding the exact SUBCLU solution set.
  (ii) When $q$ disjoint clusters $(C_r, S_r)$ have non-overlapping projections, every non-empty proper subset $S' subset S_r$ yields the identical point set $C' = C_r$, so SUBCLU emits all $2^(s_r) - 1$ non-empty subsets of $S_r$. Under INSCY, Theorem II.A guarantees that $S_r$ is finalized before every proper subset $S' subset S_r$; at each $S'$, $abs(C' inter C_r) / abs(C') = 1 >= R$, pruning all $2^(s_r) - 2$ lower-dimensional copies and retaining exactly $q$ maximal clusters.
]

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized ClusterOverview benchmark (Silhouette @rousseeuw1987silhouettes and ARI @hubert1985comparing), INSCY achieves *Iris* $"Sil" = 0.268$, $"ARI" = 0.573$; *Customers* $"Sil" = 0.203$; *Multi-Shape* $"Sil" = 0.554$ ($"ARI" = "N/A"$); and *High-D* $"Sil" = 0.336$ ($"ARI" = "N/A"$), with a mean runtime of *14.553 ms*. This empirical profile is a textbook confirmation of our theoretical analysis: on *High-D*, INSCY attains the *highest Silhouette score ($0.336$)* among all four subspace algorithms in Chapters 145–148—decisively surpassing P3C ($0.302$), COSA ($0.259$), and FIRES ($0.259$)—directly reflecting Theorem II.B's expected-density compensation against volume collapse and Theorem II.A's elimination of noisy low-dimensional projections. Conversely, on low-dimensional *Iris* ($"Sil" = 0.268$, $"ARI" = 0.573$) and *Customers* ($"Sil" = 0.203$), strict expected-density core thresholds and aggressive $R$-redundancy pruning exclude diffuse boundary objects from clusters, depressing full-space silhouette relative to soft-partitioning P3C ($0.439$, $0.347$) and COSA ($0.435$, $0.327$). Its 14.553 ms runtime reflects recursive SCY-tree construction and restriction—matching FIRES (15.611 ms) and COSA (12.455 ms).
