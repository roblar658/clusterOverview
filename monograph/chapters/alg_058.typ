#import "../template/ieee.typ": *

#alg_entry(
  id: 58,
  name: "Average Linkage Clustering (UPGMA)",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "1958",
  paradigm: "Agglomerative Unweighted Pair-Group Mean (UPGMA)",
  complexity: $O(n^2)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Introduced by Sokal and Michener @sokal1958statistical for numerical taxonomy and phylogenetic reconstruction, *Average Linkage* — formally the *Unweighted Pair-Group Method with Arithmetic Mean* (UPGMA) — resolves the extremal brittleness of Single Linkage (Alg. 56) and Complete Linkage (Alg. 57) by replacing minimum and maximum pairwise distances with their uniform cross-cluster expectation. Given a finite dissimilarity space $(X, d)$ with $X = {x_1, dots, x_n}$, the UPGMA dissimilarity between two disjoint nonempty clusters $A, B subset.eq X$ is
$
  d_U (A, B) = 1 / (abs(A) abs(B)) sum_(x in A) sum_(y in B) d(x, y).
$
The qualifier "unweighted" signifies that every *point* (or pairwise interaction $(x, y) in A times B$) receives equal weight $1 / (abs(A) abs(B))$, regardless of the internal branching topology of $A$ and $B$. Consequently, when $A$ and $B$ merge into $U = A union B$, the updated dissimilarity to any third active cluster $C subset.eq X without U$ obeys the Lance--Williams recurrence @lance1967general with size-proportional convex weights $alpha_A = abs(A) / (abs(A) + abs(B))$, $alpha_B = abs(B) / (abs(A) + abs(B))$, and $beta = gamma = 0$:
$
  d_U (A union B, C) &= abs(A) / (abs(A) + abs(B)) d_U (A, C) \
  &quad + abs(B) / (abs(A) + abs(B)) d_U (B, C).
$
This contrasts with *Weighted* Pair-Group Method with Arithmetic Mean (WPGMA / McQuitty's method), which sets fixed coefficients $alpha_A = alpha_B = 1 / 2$ irrespective of cluster cardinalities, thereby up-weighting points in small subtrees exponentially with merge depth.

In dual similarity settings where a symmetric non-negative affinity matrix $W = (w_(i j))$ replaces dissimilarities, UPGMA greedily merges the pair maximizing mean cross-similarity $s_U (A, B) = w(A, B) / (abs(A) abs(B))$, where $w(A, B) = sum_(i in A, \, j in B) w_(i j)$. Remarkably, this simple greedy rule provides a provable constant-factor global approximation to the similarity-dual Dasgupta hierarchical objective @dasgupta2016cost @moseley2017approximation.

#definition(title: "UPGMA Ultrametric and Moseley--Wang Hierarchical Objective")[
  1. The *UPGMA dendrogram* $T_U$ on $(X, d)$ iteratively merges $(A_t, B_t) in arg min_(A != B) d_U (A, B)$ at height $h_t = d_U (A_t, B_t)$ for $t = 1, dots, n - 1$, inducing the cophenetic dissimilarity $u_U (x, y) = h_(t(x, y))$ where $t(x, y)$ is the merge step uniting $x$ and $y$.
  2. Given non-negative pairwise similarities $W = (w_(i j))$, the *Moseley--Wang revenue objective* @moseley2017approximation (dual to Dasgupta's cost @dasgupta2016cost) of a rooted binary tree $T$ with leaf set $X$ is
     $
       F_"MW" (T) = sum_(1 <= i < j <= n) w_(i j) (n - abs("leaves"(T[i or j]))),
     $
     where $T[i or j]$ denotes the subtree rooted at the lowest common ancestor of $i$ and $j$. High-similarity pairs are rewarded for merging deep inside small subtrees.
]

#algbox(
  title: "Exact O(n^2) UPGMA via Nearest-Neighbor Chain",
  inputs: [Symmetric dissimilarity matrix $D in bb(R)_(>= 0)^(n times n)$ on $X = {x_1, dots, x_n}$; tie-breaking order $prec$; target $k$.],
  outputs: [Exact UPGMA merge tree ${(A_t, B_t, h_t)}_(t=1)^(n-1)$ and $k$-partition $cal(C)_k$.],
)[
  1. Initialize $cal(C) <- {{x_1}, dots, {x_n}}$, cluster sizes $n_({x_i}) <- 1$, $d_U <- D$, empty stack $S <- []$, $t <- 1$.
  2. *While* $abs(cal(C)) > 1$:
     - If $S$ is empty, push any $C in cal(C)$ onto $S$.
     - Let $A$ be the top of $S$ and compute $B = arg min_(C in cal(C) without {A}) d_U (A, C)$ (tie-breaking via $prec$).
     - *If* $B$ is on $S$ (necessarily immediately below $A$, forming an RNN pair):
       - Pop $A, B$ from $S$; record $(A_t, B_t, h_t) <- (A, B, d_U (A, B))$.
       - Form $U <- A union B$ with cardinality $n_U <- n_A + n_B$ and update $cal(C) <- (cal(C) without {A, B}) union {U}$.
       - For every $C in cal(C) without {U}$, compute
         $d_U (U, C) <- (n_A d_U (A, C) + n_B d_U (B, C)) / n_U$.
       - Increment $t <- t + 1$.
     - *Else:* push $B$ onto $S$.
  3. Cut the sorted merge sequence after $n - k$ merges to output $cal(C)_k$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

#subsec("A", "Convex Reducibility, Monotonicity, and Affine Equivariance")

#lemma(title: "Two-Sided Convex Bracket, Zero Inversions, and Affine Equivariance", origin: "Original Proof")[
  For all pairwise disjoint nonempty subsets $A, B, C subset.eq X$, UPGMA satisfies:
  1. *Convex bracket and Reducibility:*
     $
       min { d_U (A, C), \, d_U (B, C) } &<= d_U (A union B, C) \
       &<= max { d_U (A, C), \, d_U (B, C) }.
     $
  2. *Monotonicity and NN-Chain exactness:* Merge heights are non-decreasing ($h_1 <= dots <= h_(n-1)$), $u_U$ is a valid ultrametric, and NN-Chain computes the exact UPGMA tree in $Theta(n^2)$ time @mullner2011modern.
  3. *Positive affine equivariance (but nonlinear sensitivity):* For any $a > 0$ and $b >= 0$, $u_U (a d + b) = a u_U (d) + b$; however, $u_U$ is *not* invariant under strictly increasing nonlinear transforms $phi$.
]

#proof[
  *Part 1:* Writing $alpha = abs(A) / (abs(A) + abs(B)) in (0, 1)$, expanding the double sum yields
  $
    d_U (A union B, C) &= (sum_(x in A, \, z in C) d(x, z) + sum_(y in B, \, z in C) d(y, z)) / ((abs(A) + abs(B)) abs(C)) \
    &= alpha d_U (A, C) + (1 - alpha) d_U (B, C).
  $
  Because a convex combination lies between its minimum and maximum terms, the two-sided bracket holds identically.

  *Part 2:* The lower inequality $d_U (A union B, C) >= min { d_U (A, C), d_U (B, C) }$ is Bruynooghe's Reducibility property. Exactly as in Alg. 57, reducibility implies that if $(A_t, B_t)$ minimizes $d_U$ at step $t$ with height $h_t$, then every newly updated distance satisfies $d_U (A_t union B_t, C) >= min { d_U (A_t, C), d_U (B_t, C) } >= h_t$, forcing $h_(t+1) >= h_t$. Thus the dendrogram has zero inversions, $u_U in cal(U)(X)$, and RNN merges remain valid under any other merge, guaranteeing exactness and $Theta(n^2)$ complexity of NN-Chain via the $2n - 1$ push bound.

  *Part 3:* By linearity of averaging, for $d'(x, y) = a d(x, y) + b$ ($x != y$),
  $
    d'_U (A, B) &= 1 / (abs(A) abs(B)) sum_(x in A, \, y in B) (a d(x, y) + b) \
    &= a d_U (A, B) + b.
  $
  Since $a > 0$, all merge-order comparisons are preserved identically and every merge height transforms via $h'_t = a h_t + b$. To prove nonlinear non-invariance, take $X = {x_1, x_2, x_3, x_4}$ with $d(x_1, x_2) = 1$, $d(x_1, x_3) = 2$, $d(x_2, x_3) = 6$, and $d(x_1, x_4) = d(x_2, x_4) = 4.2$, $d(x_3, x_4) = 7$. After merging $A = {x_1, x_2}$ at $1$, under $d$ we have $d_U (A, {x_3}) = (2 + 6)/2 = 4 < 4.2 = d_U (A, {x_4})$, merging $x_3$ next. Under the strictly increasing transform $phi(t) = t^2$, however, $d_(U, phi) (A, {x_3}) = (4 + 36)/2 = 20 > 17.64 = d_(U, phi) (A, {x_4})$, reversing the merge order to $x_4$.
]

#subsec("B", "Exact Molecular-Clock Recovery and L_infinity Noise Robustness")

A core statistical advantage of UPGMA over extremal linkages is its exact recovery of ultrametric trees together with variance-shrinking noise resilience.

#theorem(title: "Exact Ultrametric Tree Recovery and Perturbation Stability Margin", origin: "Original Proof")[
  Let $u^(star) in cal(U)(X)$ be a true ultrametric tree with strictly positive internal level gap
  $
    Delta_"min" &= min { abs(u^(star) (x, z) - u^(star) (x, y)) : \
    &quad quad quad u^(star) (x, z) > u^(star) (x, y) } > 0.
  $
  Suppose we observe a perturbed dissimilarity $d$ with $norm(d - u^(star))_infinity = max_(x, y) abs(d(x, y) - u^(star) (x, y)) <= eta$. Then:
  1. If $eta = 0$, UPGMA recovers the exact rooted tree topology and exact ultrametric $u_U = u^(star)$ @sokal1958statistical.
  2. If $eta < Delta_"min" / 2$, UPGMA *still* recovers the exact rooted clade hierarchy of $u^(star)$ and satisfies $norm(u_U - u^(star))_infinity <= eta$.
  3. Moreover, if $d(x, y) = u^(star) (x, y) + epsilon_(x y)$ with independent zero-mean noise of variance $sigma^2$, the merge-height estimate of any true sibling clades $(A, B)$ is unbiased with variance $"Var"(d_U (A, B)) = sigma^2 / (abs(A) abs(B))$.
]

#proof[
  We prove Parts 1 and 2 simultaneously by induction on the merge step $t$. Call a subset $C subset.eq X$ a *true clade* of $u^(star)$ if it equals a maximal ultrametric ball $B_(u^(star)) (x, r)$ for some $x in X, r >= 0$. At $t = 0$, every singleton is a true clade. Assume all active clusters at step $t$ are true clades of $u^(star)$. Because $u^(star)$ is an ultrametric, for any two disjoint true clades $A, B$, $u^(star) (x, y)$ takes a single constant value $h^(star) (A, B)$ for all $(x, y) in A times B$. Averaging the pointwise error bound $abs(d(x, y) - h^(star) (A, B)) <= eta$ over $A times B$ yields
  $
    abs(d_U (A, B) - h^(star) (A, B)) &<= 1 / (abs(A) abs(B)) sum_(x in A, \, y in B) abs(d(x, y) - h^(star) (A, B)) \
    &<= eta.
  $
  Now let $(A, B)$ be any active pair attaining the minimum true merge height $h_0 = min_(P != Q) h^(star) (P, Q)$ (hence $A union B$ is a true clade of $u^(star)$). Any active pair $(P, Q)$ whose union is *not* a true clade at level $h_0$ must satisfy $h^(star) (P, Q) >= h_0 + Delta_"min"$. Consequently, when $eta < Delta_"min" / 2$,
  $
    d_U (A, B) &<= h_0 + eta < h_0 + Delta_"min" - eta \
    &<= h^(star) (P, Q) - eta <= d_U (P, Q).
  $
  Thus UPGMA strictly prefers true sibling clades at height $h_0$ over any non-clade merge, preserving the inductive invariant that every merged cluster is a true clade of $u^(star)$. When $eta = 0$, $d_U (A, B) = h^(star) (A, B)$ exactly (Part 1); when $eta < Delta_"min" / 2$, every pairwise lowest common ancestor matches that of $u^(star)$, giving $abs(u_U (x, y) - u^(star) (x, y)) = abs(d_U (A, B) - h^(star) (A, B)) <= eta$ (Part 2).

  *Part 3:* For true sibling clades $A, B$ with constant true distance $h^(star) (A, B)$,
  $
    d_U (A, B) = h^(star) (A, B) + 1 / (abs(A) abs(B)) sum_(x in A, \, y in B) epsilon_(x y).
  $
  By linearity of expectation and independence across the $abs(A) abs(B)$ distinct cross-pairs, $bb(E)[d_U (A, B)] = h^(star) (A, B)$ and $"Var"(d_U (A, B)) = sigma^2 / (abs(A) abs(B))$, showing $O(1 / (abs(A) abs(B)))$ variance attenuation at upper tree levels (in contrast to $O(sigma^2)$ extreme-value bias in Single and Complete linkage).
]

#subsec("C", "Global 1/3-Approximation for Similarity Hierarchical Clustering")

Unlike most heuristic clustering algorithms, UPGMA admits a crisp global approximation guarantee under Moseley and Wang's similarity objective @moseley2017approximation.

#theorem(title: "1/3-Approximation Guarantee of Average Linkage", origin: "Classical — Cited")[
  For any non-negative symmetric similarity matrix $W = (w_(i j))$ on $n$ items, let $T_U$ be the hierarchy returned by Average Linkage (maximizing mean similarity $s_U (A, B) = w(A, B) / (abs(A) abs(B))$) and let $T^(star) in arg max_T F_"MW" (T)$. Then @moseley2017approximation
  $
    F_"MW" (T_U) >= 1 / 3 F_"MW" (T^(star)).
  $
]

#proof[
  For any binary tree $T$ on $X$, rewrite $F_"MW" (T)$ as a sum over all unordered triples ${i, j, k} subset.eq X$. Pair ${i, j}$ contributes $w_(i j)$ to $F_"MW" (T)$ once for every leaf $k in X without "leaves"(T[i or j])$ — precisely those triples ${i, j, k}$ in which $i$ and $j$ merge *before* joining $k$. Denote the sum of similarities between $A$ and $B$ by $w(A, B) = sum_(i in A, \, j in B) w_(i j)$. At the step when $T_U$ merges clusters $A$ and $B$, let $cal(R) = cal(C) without {A, B}$ be the set of remaining active clusters. Every triple $(i, j, k) in A times B times C$ (for $C in cal(R)$) is resolved at this merge:
  - *Revenue gained:* each such triple awards $w_(i j)$ to $F_"MW" (T_U)$, totaling
    $
      "rev"(A, B) = sum_(C in cal(R)) abs(C) w(A, B).
    $
  - *Potential revenue lost:* the competing pairs $(i, k) in A times C$ and $(j, k) in B times C$ can no longer separate $j$ or $i$ first, forfeiting
    $
      "loss"(A, B) = sum_(C in cal(R)) (abs(B) w(A, C) + abs(A) w(B, C)).
    $
  Because every unordered triple ${i, j, k}$ is resolved at a unique merge of $T_U$, summing over all merges gives
  $
    sum_(A, B) ("rev"(A, B) + "loss"(A, B)) &= sum_(1 <= i < j <= n) (n - 2) w_(i j) \
    &>= F_"MW" (T^(star)),
  $
  since every pair $(i, j)$ separates at most $n - 2$ external leaves in any tree $T^(star)$. Now apply UPGMA's greedy optimality at merge $(A, B)$: for every $C in cal(R)$,
  $
    w(A, B) / (abs(A) abs(B)) >= w(A, C) / (abs(A) abs(C)), quad
    w(A, B) / (abs(A) abs(B)) >= w(B, C) / (abs(B) abs(C)).
  $
  Cross-multiplying yields $abs(C) w(A, B) >= abs(B) w(A, C)$ and $abs(C) w(A, B) >= abs(A) w(B, C)$. Adding these inequalities and summing over all $C in cal(R)$ produces
  $
    2 "rev"(A, B) >= "loss"(A, B).
  $
  Summing over all $n - 1$ merges of $T_U$ gives $F_"MW" (T^(star)) <= sum ("rev" + "loss") <= 3 sum "rev" = 3 F_"MW" (T_U)$, establishing $F_"MW" (T_U) >= (1 / 3) F_"MW" (T^(star))$.
]

#subsec("D", "Failure under Violated Molecular Clock (Long-Branch Attraction)")

While Theorem 2 establishes exact consistency on ultrametrics, UPGMA fails systematically on non-ultrametric *additive* tree metrics when branch rates vary across lineages.

#proposition(title: "Constructive Failure on Non-Ultrametric Additive Trees", origin: "Original Counterexample")[
  Consider the four-taxon additive tree $T$ on $X = {a, b, c, d}$ with true unrooted split $a b | c d$, strictly positive internal branch $e > 0$, and asymmetric leaf branch lengths $ell_a = ell_c = 1$, $ell_b = ell_d = L$. Then:
  1. The path-length metric $d_T$ satisfies the four-point condition for all $L > 0$, so $T$ is a uniquely identifiable additive phylogenetic tree.
  2. Whenever $L > 1 + e$, UPGMA erroneously merges the non-sibling short-branch pair ${a, c}$ first, recovering the false topology $a c | b d$.
]

#proof[
  The six pairwise tree distances are $d_T (a, b) = d_T (c, d) = 1 + L$, $d_T (a, c) = 2 + e$, $d_T (a, d) = d_T (b, c) = 1 + e + L$, and $d_T (b, d) = 2 L + e$. Check Buneman's four-point condition:
  $
    d_T (a, b) + d_T (c, d) &= 2 + 2 L < 2 + 2 L + 2 e \
    &= d_T (a, c) + d_T (b, d) = d_T (a, d) + d_T (b, c),
  $
  confirming $d_T$ is strictly additive with unique true split $a b | c d$. However, for $L > 1 + e$, the strictly smallest pairwise distance is $d_T (a, c) = 2 + e < 1 + L = d_T (a, b) = d_T (c, d)$ (and strictly smaller than $1 + e + L$ and $2 L + e$). Hence UPGMA merges ${a, c}$ at Step 1, destroying the true split $a b | c d$. Because UPGMA assumes equidistant root-to-leaf sums (molecular clock), slow-evolving non-sister taxa $a$ and $c$ attract; rate-corrected *Neighbor-Joining* @saitou1987neighbor resolves this failure by subtracting mean taxon divergences via the $Q$-criterion.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Computational and Memory Profile:* Using NN-Chain (Algorithm Box), exact UPGMA executes in $Theta(n^2)$ time and $Theta(n^2)$ memory (since $d_U (A, B)$ averages over all cross-pairs and cannot be updated in $Theta(n)$ space without recomputing distances). Heap-based variants run in $O(n^2 log n)$ time. For massive datasets ($n > 10^5$), exact $Theta(n^2)$ space becomes prohibitive; scalable approximate hierarchies such as Grinch @monath2019grinch and SCC @monath2021scc use nearest-neighbor graph pruning and incremental tree rotations to approximate UPGMA in near-linear time.

*Geometric Inductive Bias:* UPGMA strikes a principled balance between connectivity and compactness. Because its convex-combination update is linear in pairwise dissimilarities, it resists single-point bridge chaining (unlike SLINK) while avoiding extreme-outlier diameter jumps (unlike CLINK). It works seamlessly on *non-Euclidean* dissimilarities (Bray--Curtis, cosine, edit distance, phylogenetic substitution matrices) where Ward's centroid formula is undefined.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [*Property*], [*Characterization*],
    [Time Complexity], [$Theta(n^2)$ exact via NN-Chain @mullner2011modern; $O(n^2 log n)$ via priority queues.],
    [Space Complexity], [$Theta(n^2)$ pairwise dissimilarity matrix (no $Theta(n)$ streaming reduction known).],
    [Determinism], [Fully deterministic under fixed tie-breaking; stable under small noise (Theorem 2).],
    [Invariances], [Positive affine equivariance ($a d + b$); *not* invariant to nonlinear monotone maps.],
    [Cluster Geometry], [Moderately compact, balanced clusters; intermediate between SLINK and CLINK.],
    [Monotonicity], [Strictly inversion-free by convex reducibility; yields valid ultrametric $u_U$.],
    [Outlier Sensitivity], [Robust to isolated bridges; variance attenuates as $sigma^2 / (abs(A) abs(B))$.],
    [Primary Bottleneck], [$Theta(n^2)$ RAM footprint for large $n$ and sensitivity to rate heterogeneity.],
  ),
  caption: [Operational and theoretical profile of Average Linkage Clustering (UPGMA).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Molecular Phylogenetics under Constant-Rate Evolution:* In its original domain @sokal1958statistical, UPGMA reconstructs rooted ultrametric phylogenies directly from aligned DNA or protein distance matrices whenever mutation rates approximately satisfy a molecular clock, combining exact recovery with $O(1 / (abs(A) abs(B)))$ sampling-noise averaging (Theorem 2).
2. *Ecological Community Ordination (Bray--Curtis / Jaccard Metrics):* Community ecologists cluster species-abundance sites using non-Euclidean semi-metrics such as Bray--Curtis or UniFrac dissimilarity. Because Ward's method requires Euclidean inner products, UPGMA is the canonical variance-stabilized hierarchical clustering method in ecology and biogeography.
3. *Transcriptomic Co-Expression Heatmaps and Bioinformatics:* In bulk and single-cell RNA-seq analyses, genes are clustered using Pearson correlation distance $d = 1 - rho$. Average linkage suppresses noisy single-gene expression spikes (which break Complete linkage) while preventing unrelated pathways from chaining together (which corrupts Single linkage).
4. *Speaker Diarization and Biometric Identity Clustering:* Agglomerative Hierarchical Clustering (AHC) on PLDA log-likelihood-ratio or cosine embeddings of speech segments relies on UPGMA average similarity, directly leveraging the Moseley--Wang $1 / 3$-approximation guarantee (Theorem 3) to maximize global pairwise speaker purity.
5. *Automated Document Taxonomy and Ontology Discovery:* Building multi-level topical hierarchies from cosine similarity matrices favors UPGMA because its merge heights retain directly interpretable mean intra-cluster similarity values at every cut level.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

UPGMA's chief failure modes are (i) heterogeneous evolutionary rates violating ultrametricity (Proposition 1), (ii) strongly intertwined non-convex manifolds (where cross-cluster averages exceed intra-manifold chords), (iii) high-dimensional Euclidean concentration @beyer1999nearest (where mean pairwise distances concentrate tightly around the global mean, weakening separation contrast), and (iv) datasets with $n >> 10^5$ exceeding quadratic RAM budgets.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [*Domain*], [*Mathematical Fit*], [*Caution / Alternative*],
    [Clock-Like Phylogenetics], [Exact ultrametric recovery + noise averaging], [Rate variation; use Neighbor-Joining @saitou1987neighbor],
    [Ecological / Non-Euclidean], [Valid for arbitrary metrics without centroids], [Extreme size imbalance; check WPGMA],
    [Speaker / Embedding AHC], [$1 / 3$-optimal Moseley--Wang global revenue], [$n > 10^5$ RAM; use SCC @monath2021scc],
    [Non-Convex Manifolds], [Averages blur curved gaps across bends], [Prefer SLINK (Alg. 56) or HDBSCAN\*],
    [Euclidean Gaussian Blobs], [Solid baseline, but ignores cluster size inertia], [Prefer Ward (Alg. 59) for higher ARI],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for Average Linkage (UPGMA).],
)

#sec("V", "Structural Relations and Empirical Behavior")

Average linkage connects the non-Euclidean metric linkages (Single, Complete) to Euclidean variance minimization (Ward's Method, Alg. 59) through an exact algebraic variance decomposition.

#proposition(title: "Squared-Euclidean Decomposition and Comparison to Ward's Method", origin: "Original Derivation")[
  Let $A, B subset bb(R)^p$ be finite disjoint point sets with centroids $mu_A, mu_B$ and within-cluster sums of squares $"SSE"(A) = sum_(x in A) norm(x - mu_A)_2^2$. Then:
  1. For any dissimilarity $d$, at every step on the same active partition, $d_S (A, B) <= d_U (A, B) <= d_C (A, B)$.
  2. Under squared Euclidean dissimilarity $d^(2) (x, y) = norm(x - y)_2^2$, average linkage satisfies the exact identity
     $
       d_U^((2)) (A, B) &= norm(mu_A - mu_B)_2^2 \
       &quad + 1 / abs(A) "SSE"(A) + 1 / abs(B) "SSE"(B).
     $
  Comparing with Ward's merge increment $Delta(A, B) = (abs(A) abs(B) / (abs(A) + abs(B))) norm(mu_A - mu_B)_2^2$ (Alg. 59) reveals that squared-Euclidean UPGMA penalizes internal cluster variances additively without harmonic size weighting, whereas Ward weights squared centroid distance by harmonic cluster mass.
]

#proof[
  *Part 1:* Because the arithmetic mean of a finite set of real numbers lies between its minimum and maximum, $min_(x in A, y in B) d(x, y) <= d_U (A, B) <= max_(x in A, y in B) d(x, y)$ holds identically.

  *Part 2:* Write $x - y = (x - mu_A) + (mu_A - mu_B) - (y - mu_B)$. Expanding the squared Euclidean norm yields
  $
    norm(x - y)_2^2 &= norm(x - mu_A)_2^2 + norm(mu_A - mu_B)_2^2 + norm(y - mu_B)_2^2 \
    &quad + 2 chevron.l x - mu_A, \, mu_A - mu_B chevron.r \
    &quad - 2 chevron.l y - mu_B, \, mu_A - mu_B chevron.r - 2 chevron.l x - mu_A, \, y - mu_B chevron.r.
  $
  Averaging over all $(x, y) in A times B$, all three inner-product cross terms vanish identically because $sum_(x in A) (x - mu_A) = 0$ and $sum_(y in B) (y - mu_B) = 0$. The remaining three terms average to $"SSE"(A) / abs(A) + norm(mu_A - mu_B)_2^2 + "SSE"(B) / abs(B)$, proving the identity.
]

UPGMA also stands in close relation to WPGMA (equal-weight update $alpha_A = alpha_B = 1 / 2$), Centroid linkage (which drops the internal variance terms above and thereby loses reducibility, permitting inversions), and Neighbor-Joining @saitou1987neighbor (which corrects UPGMA's molecular-clock bias on general additive trees).

*Empirical Behavior on the ClusterOverview 4-Dataset Benchmark:* Across the standardized four-dataset suite, Average Linkage (UPGMA) achieves *Iris* $"Sil" = 0.489$, $"ARI" = 0.560$; *Customers* $"Sil" = 0.334$; *Multi-Shape* $"Sil" = 0.452$ ($"ARI" = "N/A"$); *High-D* $"Sil" = 0.245$ ($"ARI" = "N/A"$); with a mean runtime of *2.549 ms*. The empirical profile mirrors the theoretical bracket $d_S <= d_U <= d_C$ with remarkable precision. On *Multi-Shape*, UPGMA's Silhouette ($0.452$) sits squarely between Complete ($0.420$) and Single ($0.609$) linkage. On *Customers*, mean-distance averaging completely eliminates SLINK's chaining collapse ($0.067$) and outperforms Complete linkage ($0.334$ vs. $0.323$), coming within $0.013$ of Ward's best score ($0.347$). On *Iris*, UPGMA ties SLINK for the highest Silhouette ($0.489$, ahead of Complete $0.444$ and Ward $0.441$), while its $"ARI"$ ($0.560$) trails Ward ($0.696$) because equal-weight pairwise averaging does not aggressively split overlapping Gaussian petals without Ward's size-weighted inertia penalty. On *High-D*, UPGMA's Silhouette ($0.245$) is the lowest of the four linkages, consistent with the concentration of pairwise arithmetic averages around the global mean distance in high dimensions @beyer1999nearest. Finally, its mean runtime of $2.549$ ms closely matches SLINK ($2.343$ ms) and Ward ($2.344$ ms) and runs more than $2.1 times$ faster than Complete linkage ($5.391$ ms), confirming the efficiency of exact $Theta(n^2)$ NN-Chain UPGMA.
