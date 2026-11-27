#import "../template/ieee.typ": *

#alg_entry(
  id: 79,
  name: "SCC: Sub-Cluster Component Hierarchical Clustering",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2021",
  paradigm: "Round-Based Reciprocal Nearest-Neighbor Hierarchy",
  complexity: $O(N d log N + R_"rnd" |E|)$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Exact sequential hierarchical agglomerative clustering (HAC)—including the $O(N^2)$ nearest-neighbor chain algorithm @mullner2011modern—performs $N - 1$ strictly sequential merge steps, creating an insurmountable synchronization bottleneck for billion-scale distributed datasets. At the opposite extreme, parameter-free 1-NN graph methods such as FINCH @sarfraz2019finch contract all 1-NN connected components unconditionally in a single step, inducing severe chaining across weak bridges and losing monotonic linkage control. *SCC* (*Sub-Cluster Component*) @monath2021scc bridges these regimes via a *round-based, threshold-scheduled* parallel hierarchy built on reciprocal nearest neighbors (RNNs) and bounded weakly connected components of sparse similarity graphs.

Let $X = \{bold(x)_1, dots.c, bold(x)_N\} subset.eq bb(R)^d$ and let $G_0 = (V_0, E_0, w_0)$ be a sparse undirected similarity graph on $V_0 = X$ (e.g., a symmetric $k_"nn"$ graph with cosine or RBF affinities $w_0 (bold(x), bold(y)) > 0$). SCC takes a monotonically decreasing sequence of $R$ linkage similarity thresholds
$ tau_1 > tau_2 > dots.c > tau_R >= 0, $
together with a symmetric inter-cluster linkage function $f(A, B)$ (typically average linkage $f_"avg" (A, B) = (|A| |B|)^(-1) sum_(x in A, y in B) w_0 (x, y)$, or single/complete linkage). At round $r in \{1, dots.c, R\}$, SCC maintains a coarsened cluster graph $G_r = (V_r, E_r, f)$ where each vertex $u in V_r$ represents a disjoint subset $C_u subset.eq X$ and each edge $\{u, v\} in E_r$ carries exact inter-cluster similarity $f(C_u, C_v)$.

#definition(title: "Thresholded 1-NN Digraph, RNN Matching, and SCC Contraction")[
  Given cluster graph $G_r = (V_r, E_r, f)$ with a total tie-breaking order on vertex IDs $"id": V_r -> bb(N)$ and threshold $tau_r$:
  1. For each vertex $u in V_r$, its *thresholded best neighbor* $"nn"_(tau_r) (u)$ is defined iff $max_(v in cal(N)_r (u)) f(u, v) >= tau_r$ as the unique vertex $v^* in cal(N)_r (u)$ maximizing the lexicographic pair $(f(u, v), -"id"(v))$.
  2. The *thresholded 1-NN digraph* is $H_(tau_r) = (V_r, A_(tau_r))$ with directed arc set
     $ A_(tau_r) = \{ (u, "nn"_(tau_r) (u)) : &u in V_r, \
       &max_(v in cal(N)_r (u)) f(u, v) >= tau_r \}. $
  3. An unordered pair $\{a, b\} subset.eq V_r$ is a *Reciprocal Nearest-Neighbor (RNN) pair* at threshold $tau_r$ iff both $(a, b) in A_(tau_r)$ and $(b, a) in A_(tau_r)$. Denote the set of all RNN pairs by $cal(M)_(tau_r)$.
  4. In *strict RNN mode*, SCC contracts every pair in $cal(M)_(tau_r)$ simultaneously; in *sub-cluster component (SCC) mode*, SCC contracts every weakly connected component $K$ of $H_(tau_r)$ simultaneously (optionally capping component merge diameter or iterating sub-rounds until $A_(tau_r) = emptyset$).
]

After each contraction step, edge weights between newly contracted super-nodes $U, W in V_(r+1)$ are updated in parallel using exact additive sufficient statistics—specifically $W(U, W) = sum_(u in U, v in W) W(u, v)$ and $f_"avg" (U, W) = W(U, W) / (|C_U| |C_W|)$ (or exact vector sums $(n_U, bold(s)_U)$ when $w_0$ is an inner product).

#algbox(
  title: "SCC (Scalable Sub-Cluster Component Hierarchical Clustering)",
  inputs: [Sparse similarity graph $G_0 = (V_0, E_0, w_0)$, threshold schedule $tau_1 > dots.c > tau_R$, linkage $f$, contraction mode (RNN or Full-Component).],
  outputs: [Multi-level hierarchy $(cal(P)_0, cal(P)_1, dots.c, cal(P)_R)$ of nested partitions of $X$.],
)[
  1. Initialize $V_1 <- V_0$ (singleton clusters), $E_1 <- E_0$, and partition $cal(P)_0 <- \{\{bold(x)\} : bold(x) in X\}$.
  2. *For* round $r = 1, dots.c, R$:
     - *Repeat* (sub-rounds at fixed threshold $tau_r$):
       - *Parallel Best-Neighbor Step*: For every $u in V_r$ in parallel, compute $"nn"_(tau_r) (u)$ among incident edges with $f(u, v) >= tau_r$ (breaking ties by minimum $"id"(v)$), forming digraph $H_(tau_r) = (V_r, A_(tau_r))$.
       - *If* $A_(tau_r) = emptyset$: *break* (threshold $tau_r$ exhausted).
       - *Parallel Contraction*:
         - *If* strict RNN mode: extract all 2-cycles $cal(M)_(tau_r) = \{\{a, b\} : (a, b), (b, a) in A_(tau_r)\}$ and contract each pair into a super-node with ID $max("id"(a), "id"(b))$.
         - *If* full SCC mode: compute weakly connected components $\{K_1, dots.c, K_m\}$ of $H_(tau_r)$ via parallel union-find / label propagation and contract each $K_j$ into a single super-node.
       - *Parallel Graph Re-Weighting*: Aggregate inter-component edges in parallel using exact linkage update $f(U, W)$ (maintaining $W(U, W)$ and $|C_U|$ for UPGMA, or $min / max$ for complete/single linkage).
     - Record level-$r$ flat partition $cal(P)_r$ and set $G_(r+1) <- G_r$.
  3. *Return* nested hierarchy $(cal(P)_0, cal(P)_1, dots.c, cal(P)_R)$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five rigorous theoretical results governing SCC: a complete graph-theoretic decomposition of thresholded 1-NN digraphs into RNN-rooted pseudotrees (valid even under arbitrary weight ties), order-invariance of parallel RNN contraction relative to exact sequential HAC under reducible linkages, worst-case distortion and logarithmic round-depth bounds for full-component contraction, exact cluster recovery under thresholded separability, and a constructive 1D counterexample characterizing chaining failure under overly coarse schedules.

#subsec("A", "Unique RNN Decomposition of Thresholded 1-NN Digraphs")

#theorem(title: "Unique 2-Cycle Pseudotree Structure Under Tie-Breaking", origin: "Original Proof")[
  Let $G = (V, E, f)$ be any finite undirected weighted graph with injective vertex IDs $"id": V -> bb(N)$, and let $H_tau = (V, A_tau)$ be the thresholded 1-NN digraph from Definition I.1. Then:
  1. Every vertex $v in V$ with in-degree $d^-(v) >= 1$ has out-degree $d^+(v) = 1$.
  2. Every non-trivial weakly connected component $K$ of $H_tau$ ($|K| >= 2$) contains *exactly one* directed cycle, and that cycle has length exactly $2$—a unique Reciprocal Nearest-Neighbor pair $\{a, b\} subset.eq K$ achieving the maximum arc weight in $K$.
  3. Deleting the two arcs $(a, b)$ and $(b, a)$ partitions $K$ into two vertex-disjoint directed in-trees rooted at $a$ and $b$, along whose paths edge similarities are monotonically non-decreasing.
]
#proof[
  *Part 1:* Suppose $d^-(v) >= 1$, so there exists $(u, v) in A_tau$. By Definition I.1, $f(u, v) >= tau$. By symmetry of $f$, $v$ has at least one neighbor ($u$) with $f(v, u) = f(u, v) >= tau$, hence $max_(w in cal(N)(v)) f(v, w) >= tau$ and $v$ emits a unique outgoing arc ($d^+(v) = 1$).

  *Part 2:* Let $K$ be any weakly connected component with $|K| >= 2$. Every vertex $v in K$ has total degree $d^-(v) + d^+(v) >= 1$: if $d^+(v) = 0$, then $d^-(v) >= 1$, contradicting Part 1. Hence *every* $v in K$ satisfies $d^+(v) = 1$, making $H_tau [K]$ a connected functional digraph with $|K|$ vertices and $|K|$ arcs. Every connected functional digraph contains a unique directed cycle $C = (u_0 -> u_1 -> dots.c -> u_(L-1) -> u_0)$ with $L >= 2$.

  We prove $L = 2$ even when arbitrary edge weights are tied. Because each $u_i$ preferred $u_(i+1)$ over its undirected neighbor $u_(i-1)$ (indices mod $L$), we have $f(u_i, u_(i+1)) >= f(u_i, u_(i-1))$ for all $i$. Chaining around the cycle gives
  $ f(u_0, u_1) <= f(u_1, u_2) <= dots.c <= f(u_(L-1), u_0) <= f(u_0, u_1), $
  forcing all cycle arcs to share a single weight $c$. Suppose for contradiction that $L >= 3$. Let $u_m$ be the unique vertex on $C$ minimizing $"id"(u_m)$, and consider $u_(m+1)$. Both $u_m$ and $u_(m+2)$ are neighbors of $u_(m+1)$ with equal weight $c$, and since $L >= 3$, $u_(m+2) != u_m$. By minimality of $"id"(u_m)$, $"id"(u_m) < "id"(u_(m+2))$, contradicting that $u_(m+1)$ selected $u_(m+2)$ under the lexicographic rule $(f, -"id")$. Thus $L = 2$, yielding a unique RNN pair $\{a, b\}$.

  *Part 3:* Removing the unique cycle edges $(a, b)$ and $(b, a)$ from the connected functional graph leaves $|K| - 2$ edges and zero cycles on $|K|$ vertices, i.e., a forest of two directed in-trees rooted at $a$ and $b$. For any consecutive arcs $u -> v -> w$ along an in-tree path, $v$'s choice of $w$ over $u$ gives $f(v, w) >= f(u, v)$, proving monotonic non-decreasing similarity toward $\{a, b\}$ and maximality of $f(a, b)$ among all arcs in $K$.
]

#subsec("B", "Exact Equivalence of Parallel RNN Contraction to Sequential HAC")

#theorem(title: "Order-Invariance of Parallel RNN Contraction", origin: "Classical — Cited")[
  Suppose linkage $f$ satisfies *Bruynooghe's reducibility condition* @mullner2011modern @monath2021scc: for all pairwise disjoint clusters $A, B, C$,
  $ f(A union B, C) <= max(f(A, C), f(B, C)). $
  (This holds for single, complete, average/UPGMA, and negated Ward linkages.) Then the RNN matching $cal(M)_tau$ in $H_tau$ consists of pairwise vertex-disjoint pairs, and contracting all pairs in $cal(M)_tau$ simultaneously—assigning each merged node $a union b$ the ID $max("id"(a), "id"(b))$—preserves the RNN property of every remaining pair and produces the exact same hierarchy (up to equal-similarity tie-breaks) as exact sequential HAC down to threshold $tau$.
]
#proof[
  By Theorem II-A, every vertex belongs to at most one 2-cycle in $H_tau$, so $cal(M)_tau = \{\{a_1, b_1\}, dots.c, \{a_m, b_m\}\}$ is a vertex-disjoint matching. Consider contracting $\{a_1, b_1\}$ into $M_1 = a_1 union b_1$, and fix any other pair $\{a_j, b_j\} in cal(M)_tau$ ($j >= 2$). By reducibility and the fact that $b_j = "nn"_tau (a_j)$,
  $ f(M_1, a_j) &= f(a_1 union b_1, a_j) \
                &<= max(f(a_1, a_j), f(b_1, a_j)) <= f(a_j, b_j). $
  If $f(M_1, a_j) < f(a_j, b_j)$, $a_j$ strictly prefers $b_j$ over $M_1$. If equality holds, then without loss of generality $f(a_1, a_j) = f(a_j, b_j)$; since $a_j$ originally preferred $b_j$ over $a_1$, we had $"id"(b_j) < "id"(a_1) <= max("id"(a_1), "id"(b_1)) = "id"(M_1)$, so $a_j$ *still* strictly prefers $b_j$ over $M_1$! The identical argument applies to $b_j$, proving $\{a_j, b_j\}$ remains an RNN pair after contracting $\{a_1, b_1\}$. Inducting sequentially over $i = 1, dots.c, m$ shows that parallel contraction of $cal(M)_tau$ executes the exact same merge set as $m$ valid nearest-neighbor-chain merges in sequential HAC @mullner2011modern.
]

#subsec("C", "Approximation Ratio and Logarithmic Round Complexity of Full SCC")

#theorem(title: "Component-Contraction Bounds and Logarithmic Round Depth", origin: "Original Derivation")[
  Consider full-component SCC executing with a geometric schedule $tau_r = alpha tau_(r-1)$ ($alpha in (0, 1)$), exhausting each threshold $tau_(r-1)$ before entering round $r$. Then:
  1. Every contracted component $K$ at round $r$ is connected in the threshold graph $G_(>= tau_r)$ and satisfies single-linkage cohesion $min_(emptyset != S subset.neq K) max_(u in S, v in K without S) f(u, v) >= tau_r$.
  2. Every arc $(u, v) in A_(tau_r)$ contracted during round $r$ satisfies $tau_r <= f(u, v) < tau_(r-1)$ (under reducible linkage), guaranteeing a per-merge similarity ratio of at most $1 / alpha = tau_(r-1) / tau_r$ relative to the global maximum active edge.
  3. In any sub-round where a fraction $rho in (0, 1]$ of vertices have non-empty thresholded neighborhoods, the active vertex count shrinks by $|V'| <= (1 - rho / 2) |V|$; when $rho = 1$, at most $ceil(log_2 N)$ sub-rounds suffice to contract the entire graph to a single root.
]
#proof[
  *Part 1:* Every undirected edge underlying an arc of $H_(tau_r)$ has weight $>= tau_r$. Because $K$ is weakly connected in $H_(tau_r)$, any non-trivial bipartition $(S, K without S)$ of $K$ is crossed by at least one arc of $H_(tau_r)$ of weight $>= tau_r$, establishing the cut bound.

  *Part 2:* At the end of round $r - 1$, $A_(tau_(r-1)) = emptyset$, so every active edge $\{u, v\} in E_r$ satisfies $f(u, v) < tau_(r-1)$. By reducibility (applied inductively to any component contraction $K_1, K_2$), merging clusters whose pairwise similarities are all strictly below $tau_(r-1)$ can never create a new inter-cluster similarity $>= tau_(r-1)$:
  $ f(K_1, K_2) <= max_(u in K_1, v in K_2) f(u, v) < tau_(r-1). $
  Hence throughout round $r$, every contracted arc satisfies $tau_r <= f(u, v) < tau_(r-1)$, bounding the ratio between the maximum remaining edge and any contracted edge by $tau_(r-1) / tau_r = 1 / alpha$.

  *Part 3:* Let $V_"act" subset.eq V$ be the set of vertices incident to at least one arc of $H_(tau_r)$, with $|V_"act"| >= rho |V|$. By Theorem II-A, every weakly connected component inside $V_"act"$ has size $|K_j| >= 2$, so contracting $V_"act"$ yields at most $|V_"act"| / 2$ super-nodes. Including the $|V| - |V_"act"|$ isolated vertices gives
  $ |V'| &<= (|V| - |V_"act"|) + |V_"act"| / 2 \
         &= |V| - |V_"act"| / 2 <= (1 - rho / 2) |V|. $
  For $rho = 1$, $|V'| <= floor(|V| / 2)$, yielding total contraction depth at most $ceil(log_2 N)$.
]

#subsec("D", "Exact Cluster Recovery Under Thresholded Separability")

#theorem(title: "Exact Recovery Under Cluster Separability", origin: "Classical — Cited")[
  Let $cal(C)^* = \{C_1^*, dots.c, C_K^*\}$ satisfy *thresholded cluster separability* @monath2021scc: there exist $tau_"in" > tau_"out" >= 0$ such that for every $k in [K]$, every nonempty disjoint $A, B subset.eq C_k^*$, and every nonempty $S subset.eq X without C_k^*$,
  $ f(A, B) >= tau_"in" > tau_"out" >= f(A, S). $
  If the SCC schedule includes at least one threshold $tau_(r^*) in (tau_"out", tau_"in"]$ (with all intra-cluster edges retained in $G_0$), then upon completion of round $r^*$, the level partition $cal(P)_(r^*)$ equals $cal(C)^*$ identically and the resulting hierarchy satisfies $"DP"(cal(T), cal(C)^*) = 1$.
]
#proof[
  We prove by induction across all sub-rounds up to the end of round $r^*$ that every active cluster $u in V_r$ is *pure* ($C_u subset.eq C_k^*$ for some $k$). At $r = 0$, every singleton is trivially pure. Suppose every active cluster is pure prior to a sub-round at threshold $tau_r >= tau_(r^*) > tau_"out"$. For any two clusters $u, v$ belonging to distinct ground-truth classes $C_k^* != C_j^*$, separability gives $f(C_u, C_v) <= tau_"out" < tau_r$. Thus *no* cross-cluster edge can enter $H_(tau_r)$, so every weakly connected component (and every RNN pair) consists solely of clusters from a single class $C_k^*$, preserving purity after contraction.

  Now consider the end of round $r^*$, where $tau_(r^*) <= tau_"in"$ and $A_(tau_(r^*)) = emptyset$. Suppose for contradiction that some ground-truth cluster $C_k^*$ remained split across two or more pure active clusters $u, v subset.eq C_k^*$. By separability, $f(C_u, C_v) >= tau_"in" >= tau_(r^*)$, so $u$ and $v$ have an edge exceeding $tau_(r^*)$, contradicting $A_(tau_(r^*)) = emptyset$. Hence $cal(P)_(r^*) = cal(C)^*$, so every intra-cluster pair has its LCA inside its pure cluster node, yielding $"DP"(cal(T), cal(C)^*) = 1$.
]

#subsec("E", "Constructive Counterexample: Coarse-Threshold Chaining")

While Theorem II-C shows geometric schedules preserve a $1 / alpha$ approximation factor, an overly coarse single threshold $tau_1$ degrades full-component SCC via directed path chaining—even when exact average-linkage HAC and finely scheduled SCC recover the ground truth perfectly.

#proposition(title: "Component-Chaining Failure Under Coarse Threshold Schedules", origin: "Original Counterexample")[
  There exists a 4-point 1D dataset $X = \{x_1, x_2, x_3, x_4\} subset.eq bb(R)$ partitioned into $C_1^* = \{x_1, x_2\}$, $C_2^* = \{x_3, x_4\}$ under linear average-linkage similarity $f(x, y) = 2 - |x - y|$ such that:
  1. Single-threshold full-component SCC with $tau_1 = 1.25$ chains all four points into a single flat component with $"DP" = 0.50$ (and at most $0.75$ under any internal binarization preserving $\{x_1, x_2\}$).
  2. Two-threshold SCC with $tau_1 = 1.50 > tau_2 = 1.25$ (as well as strict RNN mode and exact UPGMA) recovers $\{C_1^*, C_2^*\}$ with $"DP" = 1.0$.
]
#proof[
  Set $x_1 = 0.0$, $x_2 = 0.4$, $x_3 = 1.0$, and $x_4 = 1.7$, with pairwise similarities $f(x_1, x_2) = 1.6$, $f(x_2, x_3) = 1.4$, $f(x_3, x_4) = 1.3$, and $f(x_1, x_3) = 1.0$ (all remaining pairs $<= 0.7$). The initial 1-NN arcs are $x_1 -> x_2$ ($1.6$), $x_2 -> x_1$ ($1.6$), $x_3 -> x_2$ ($1.4 > 1.3$), and $x_4 -> x_3$ ($1.3$).
  - *Coarse Schedule ($tau_1 = 1.25$)*: Every 1-NN arc exceeds $1.25$, producing a single weakly connected path $x_4 -> x_3 -> x_2 <-> x_1$. Full-component SCC contracts all four points simultaneously into one node of size $4$, yielding $"DP" = 1 / 2 (2 / 4 + 2 / 4) = 0.50$.
  - *Fine Schedule ($tau_1 = 1.50, tau_2 = 1.25$)*: At $tau_1 = 1.50$, only $x_1 <-> x_2$ ($1.6 >= 1.50$) enters $H_(tau_1)$, contracting into $M = \{x_1, x_2\}$. At $tau_2 = 1.25$, exact average-linkage re-weighting computes
    $ f_"avg" (M, \{x_3\}) &= (f(x_1, x_3) + f(x_2, x_3)) / 2 \
                           &= (1.0 + 1.4) / 2 = 1.20 < 1.25, $
    whereas $f_"avg" (\{x_3\}, \{x_4\}) = 1.30 >= 1.25$. Thus $x_3 <-> x_4$ contracts cleanly into $C_2^* = \{x_3, x_4\}$, separating $C_1^*$ and $C_2^*$ with $"DP" = 1.0$.
]

#sec("III", "Complexity, Geometry and Operational Profile")

*Time, Space, and Distributed Scalability.* Constructing an approximate $k_"nn"$ similarity graph $G_0$ with $|E_0| = O(N k_"nn")$ edges takes $O(N d log N)$ using modern ANN indices. Each SCC sub-round performs one linear pass over active edges ($O(|E_r|)$ work) to compute best neighbors, parallel connected components / RNN matching, and edge aggregation—all embarrassingly data-parallel in MapReduce, Spark, or shared-memory PRAM models @monath2021scc. With $R_"rnd"$ total sub-rounds (typically $O(R + log N)$ via Theorem II-C), total work is $O(N d log N + R_"rnd" |E_0|)$ and space is $O(N d + |E_0|)$, scaling empirically to hundreds of millions of embeddings across distributed workers.

*Geometric Invariances.* SCC operates purely on the pairwise similarity graph $G_0$ and linkage $f$. Consequently, when $G_0$ is constructed from Euclidean distances or cosine similarities, SCC inherits full *permutation invariance* (zero arrival-order bias, unlike PERCH and GRINCH!), *orthogonal rotation invariance*, and *translation invariance* (for metric affinities). Moreover, monotonic strictly increasing rescalings of affinities $g compose f$ yield identical hierarchies provided the schedule transforms identically as $g(tau_r)$.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(N d log N)$ graph build + $O(R_"rnd" |E_0|)$ parallel contraction rounds],
    [Space Complexity], [$O(N d + N k_"nn")$ sparse similarity graph and cluster map],
    [Parallelism Model], [Data-parallel bulk-synchronous rounds ($O(log N)$ depth in Theorem II-C)],
    [Invariances], [Arrival-order / permutation invariant; orthogonal invariant under cosine/$L_2$],
    [Primary Hyperparameters], [Threshold schedule $(tau_1, dots.c, tau_R)$, sparse graph degree $k_"nn"$, linkage $f$],
    [Hierarchy Format], [Multi-level $R$-layer nested hierarchy (convertible to binary dendrogram)],
  ),
  caption: [Operational and geometric profile of SCC (Alg. 79).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

1. *Billion-Scale Web Document and Entity Clustering.* Industrial search and knowledge-graph pipelines cluster $N > 10^8$ dense embeddings using distributed MapReduce jobs @monath2021scc. SCC's $O(log N)$ parallel round depth and sparse graph contraction eliminate the $O(N)$ sequential barrier of exact HAC and incremental trees.
2. *Large-Scale Social and Citation Graph Community Hierarchies.* Real-world interaction networks are already sparse graphs ($|E_0| = O(N)$); SCC directly ingests $G_0$ without vectorization, producing multi-resolution community trees aligned with natural reciprocal-affinity cores.
3. *Massive-Scale Protein Sequence and Homology Family Clustering.* Metagenomic protein databases contain hundreds of millions of sequences connected by sparse sequence-alignment scores (MMseqs2/BLAST graphs); thresholded RNN contraction groups paralogs and orthologs cleanly across evolutionary distance scales without single-linkage domain-fusion chaining.
4. *Web-Scale Near-Duplicate Image and Video Clustering.* Content-moderation and media-deduplication platforms construct sparse ANN cosine graphs over billions of frames; SCC's fine geometric threshold schedule isolates exact duplicates at high $tau_r$ and semantic event clusters at lower $tau_r$.
5. *Distributed Customer Segmentation on Sparse Interaction Graphs.* Retail and behavioral datasets with heterogeneous cluster sizes benefit from order-invariant reciprocal matching, which avoids both early-insertion lock-in and forced spherical Voronoi splits.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

SCC exhibits four principal failure modes: (i) *Overly coarse threshold spacing*, which triggers directed-path chaining in full-component mode as proven in Proposition II-E; (ii) *Sparse graph disconnection* (too small $k_"nn"$), where long-range true inter-cluster edges are pruned from $G_0$, halting merges prematurely at lower thresholds unless missing similarities are computed on-the-fly from stored vector sums; (iii) *Small-dataset orchestration overhead*, where graph construction and multi-round threshold sweeps carry higher constant latency than lightweight incremental trees like GRINCH; and (iv) *Non-reducible linkages* (centroid / median linkage), where Theorem II-B's exact parallel equivalence breaks down due to linkage inversions.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Billion-scale distributed embedding clustering], [Data-parallel $O(log N)$ rounds; sparse graph scaling], [Use fine geometric schedule ($alpha >= 0.85$)],
    [Sparse similarity / biological homology graphs], [Native graph input; avoids dense $O(N^2)$ matrix], [Increase $k_"nn"$ or fallback to vector centroids],
    [Order-sensitive scientific clustering], [100% arrival-order invariant; exact RNN equivalence], [Ensure reducible linkage (UPGMA/complete/single)],
    [Low-latency single-point streaming], [Weak fit: batch round-based graph contraction], [Use GRINCH (Alg. 77) or BICO (Alg. 78)],
    [Tiny-$N$ real-time embedded workloads], [High per-round overhead relative to $N$], [Use GRINCH (Alg. 77) or exact NN-chain HAC],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for SCC.],
)

#sec("V", "Structural Relations and Empirical Behavior")

SCC unifies exact nearest-neighbor-chain HAC @mullner2011modern, parameter-free 1-NN contraction (FINCH @sarfraz2019finch), and incremental hierarchical clustering (PERCH @kobren2017perch, Alg. 76; GRINCH @monath2019grinch, Alg. 77). Setting $R = 1$ with $tau_1 = 0$ in full-component mode recovers a single pass of FINCH's first-neighbor connected components, whereas letting the threshold spacing approach zero in strict RNN mode recovers exact sequential HAC (Theorem II-B).

#proposition(title: "Exact Reduction Between FINCH, SCC, and Sequential UPGMA", origin: "Original Proof")[
  Let $G_0$ be the complete Euclidean/cosine similarity graph on $X$ with reducible linkage $f$ and distinct pairwise weights. Then:
  1. Full-component SCC with a single threshold $tau_1 = -infinity$ (one contraction per round without sub-round re-checking) coincides identically with FINCH @sarfraz2019finch at every partition level.
  2. Strict-RNN SCC with any schedule whose distinct thresholds separate all $N - 1$ merge levels (or with a single threshold $tau_1 = -infinity$ iterating RNN sub-rounds) produces the *exact* dendrogram of sequential HAC @mullner2011modern.
]
#proof[
  For part 1, at $tau_1 = -infinity$, every vertex $u in V_r$ emits an arc to its global 1-NN, so $H_(-infinity)$ is precisely the directed 1-NN graph whose undirected weakly connected components define FINCH's adjacency contraction $A(i, j) = 1$ iff $j = kappa_i^1$ or $i = kappa_j^1$ (ignoring shared-neighbor extensions, which form a superset of the same 1-NN components). For part 2, when $tau_1 = -infinity$, Theorem II-A guarantees $cal(M)_(-infinity) != emptyset$ whenever $|V_r| >= 2$. By Theorem II-B, each parallel RNN contraction executes a valid batch of nearest-neighbor-chain merges without altering the merge order of any unaffected pair, terminating after at most $N - 1$ sub-rounds in the unique exact HAC dendrogram @mullner2011modern.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark @rousseeuw1987silhouettes @hubert1985comparing, SCC records: *Iris* Silhouette $= 0.489$, ARI $= 0.560$; *Customers* Silhouette $= 0.381$; *Multi-Shape* Silhouette $= 0.609$ (ARI N/A); *High-D* Silhouette $= 0.229$ (ARI N/A); with mean runtime $= 1.590$ ms. Notably, SCC achieves the *highest Customers Silhouette* ($0.381$) across all four algorithms (Algs. 76--79), outperforming BICO ($0.354$), GRINCH ($0.311$), and PERCH ($0.273$), because arrival-order-invariant reciprocal contraction cleanly separates heterogeneous spending profiles without spherical centroid bias or early-insertion splits. On Iris ($0.489$ / $0.560$) and Multi-Shape ($0.609$), SCC matches GRINCH to three decimal places—empirically confirming Proposition V.1 and GRINCH's convergence toward exact average-linkage HAC—while its $1.590$ ms small-dataset runtime reflects multi-round graph overhead that amortizes dramatically only at massive distributed scale.
