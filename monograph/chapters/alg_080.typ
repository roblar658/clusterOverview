#import "../template/ieee.typ": *

#alg_entry(
  id: 80,
  name: "FINCH: First Integer Neighbor Clustering Hierarchy",
  category: "Cat. III: Hierarchical & Tree-Based",
  year: "2019",
  paradigm: "Parameter-Free First-Neighbor Hierarchical Contraction",
  complexity: $O(n log n) " (ANN), " O(n^2 d) " (exact)"$,
  provenance: "Original Proofs + Cited Classical Bounds",
)

#sec("I", "Problem Formulation and Mathematical Model")

Most clustering methods require sensitive hyperparameters---such as cluster count $k$, neighborhood radius $epsilon$, or linkage cut thresholds---that strongly affect partition quality @jain2010data, @mullner2011modern. The *First Integer Neighbor Clustering Hierarchy* (FINCH) algorithm @sarfraz2019finch eliminates all continuous tuning parameters by observing that every observation's single nearest neighbor defines a parameter-free adjacency graph whose connected components immediately form a valid clustering. Iterating this contraction on cluster centroids builds a compact hierarchy from raw points to a single root.

Let $X = {x_1, dots, x_n} subset.eq bb(R)^d$ ($n >= 2$) be a dataset with metric or symmetric dissimilarity $rho : bb(R)^d times bb(R)^d -> [0, infinity)$, typically Euclidean distance $rho_2 (x, y) = norm(x - y)_2$ or cosine distance $rho_"cos" (x, y) = 1 - chevron.l x, y chevron.r / (norm(x)_2 norm(y)_2)$ on $bb(R)^d without {0}$, indexed by $V = {1, dots, n}$.

#definition(title: "First-Neighbor Adjacency and FINCH Hierarchy")[
  For each $i in V$, its *first nearest neighbor* is
  $ kappa_i^1 = arg min_(j in V without {i}) rho(x_i, x_j), $
  with deterministic tie-breaking (smallest index) when non-unique. Let $P in {0, 1}^(n times n)$ be the directed 1-NN indicator matrix with $P_(i j) = bb(I)[j = kappa_i^1]$. The symmetric *FINCH adjacency matrix* $A in {0, 1}^(n times n)$ is
  $ A(i, j) = cases(
      1 & "if" j = kappa_i^1 "or" i = kappa_j^1,
      1 & "if" kappa_i^1 = kappa_j^1,
      0 & "otherwise,"
    ) $
  equivalently $A = bb(I)[P + P^top + P P^top > 0]$. Let $cal(P)_1 = {C_(1, 1), dots, C_(n_1, 1)}$ be the connected components of $G_A = (V, E_A)$. Initializing $n_0 = n$, $mu_(i, 0) = x_i$, and $w_(i, 0) = 1$, FINCH builds nested partitions $cal(P)_1, dots, cal(P)_L$ recursively: at level $ell >= 1$ ($n_(ell - 1) >= 2$), it finds first neighbors $kappa_a^(1, ell)$ among centroids ${mu_(a, ell - 1)}$, extracts components ${C_(r, ell)}_(r = 1)^(n_ell)$, and computes weighted centroids
  $ w_(r, ell) &= sum_(a in C_(r, ell)) w_(a, ell - 1), \
    mu_(r, ell) &= 1 / w_(r, ell) sum_(a in C_(r, ell)) w_(a, ell - 1) mu_(a, ell - 1), $
  stopping at the first $L$ with $n_L = 1$ (returning $cal(P)_1, dots, cal(P)_(L - 1)$).
]

Under $rho_2$, $mu_(r, ell)$ equals the exact sample mean of all original points in cluster $r$ at level $ell$; under $rho_"cos"$, centroids may be $ell_2$-normalized before subsequent 1-NN queries @sarfraz2019finch. If a fixed target $K in.not {n_ell}$ is required, Sarfraz et al. @sarfraz2019finch select the coarsest level $ell^*$ with $n_(ell^*) >= K$ and greedily merge the closest centroid pairs down to $K$.

#algbox(
  title: "FINCH (First Integer Neighbor Clustering Hierarchy)",
  inputs: [Data $X = {x_1, dots, x_n} subset.eq bb(R)^d$, metric $rho$, optional $K$.],
  outputs: [Nested hierarchy ${cal(P)_1, dots, cal(P)_(L - 1)}$ (and optional $K$-cut).],
)[
  + Set $ell <- 0$, $n_0 <- n$, $mu_(i, 0) <- x_i$, $w_(i, 0) <- 1$, $y_i^((0)) <- i$.
  + *While* $n_ell > 1$:
    - Increment $ell <- ell + 1$ and compute $kappa_a^1 = arg min_(b != a) rho(mu_(a, ell - 1), mu_(b, ell - 1))$ via exact GEMM or tree / ANN index.
    - Find components ${C_(r, ell)}_(r = 1)^(n_ell)$ of $G_1^((ell))$ with edges ${(a, kappa_a^1)}$ via union-find (identical to $G_A^((ell))$ by Theorem 1).
    - Update labels $y_i^((ell)) <- r$ for $y_i^((ell - 1)) in C_(r, ell)$ and aggregate weighted centroids $(mu_(r, ell), w_(r, ell))$.
  + Return ${cal(P)_1, dots, cal(P)_(L - 1)}$. If target $K in.not {n_ell}$ is requested, merge closest centroid pairs from $cal(P)_(ell^*)$ ($ell^* = max {ell : n_ell >= K}$) down to $K$.
]

#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")

We establish five theoretical results for FINCH: an algebraic graph reduction and exact trace cluster count, tight logarithmic depth, the stochastic-geometry compression ratio, invariances with an asynchronous collapse counterexample, and amortized complexity with high-dimensional concentration limits.

#subsec("A", "Component Equivalence, MNN Bijection, and MST Subforest")

Although @sarfraz2019finch defines $A = bb(I)[P + P^top + P P^top > 0]$, the quadratic shared-neighbor term $P P^top$ leaves connected components unchanged while risking quadratic edge blow-up.

#theorem(title: "Component Equivalence, Trace Count, and Boruvka Subforest", origin: "Original Proof")[
  Let $(X, rho)$ be a finite metric space ($n >= 2$) with distinct pairwise distances, $G_1 = (V, E_1)$ the undirected 1-NN graph ($E_1 = {{i, kappa_i^1} : i in V}$), and $G_A = (V, E_A)$ the graph of $A = bb(I)[P + P^top + P P^top > 0]$. Then:
  + $G_A$ and $G_1$ have identical connected components, whereas $|E_1| <= n - 1$ while $|E_A|$ can reach $Theta(n^2)$.
  + Every component $C$ of $G_1$ is a tree containing a unique mutual nearest-neighbor (MNN) pair ${u, v}$ ($kappa_u^1 = v, kappa_v^1 = u$), with strictly decreasing edge weights along directed paths into ${u, v}$.
  + The level-1 cluster count satisfies the exact trace identity
    $ n_1 = 1 / 2 "Tr"(P^2) = n - |E_1|. $
  + $G_1$ is a spanning subforest of the unique MST of $(X, rho)$, coinciding with one parallel contraction round of Boruvka's MST algorithm.
]

#proof[
  *(1)* Clearly $E_1 subset.eq E_A$. If ${i, j} in E_A without E_1$, then $(P P^top)_(i j) = 1$, so $kappa_i^1 = k = kappa_j^1$ for some $k in V$. Thus ${i, k}, {j, k} in E_1$, giving a path $i - k - j$ in $G_1$. By induction on path length, $G_A$ and $G_1$ share identical connected components. For the edge blow-up, set $x_1 = 0$ and $x_i = (1 - epsilon_i) e_(i - 1) in bb(R)^(n - 1)$ ($2 <= i <= n$) with $0 < epsilon_2 < dots < epsilon_n < 1/4$. Then $rho_2 (x_i, x_1) < 1 < 3 sqrt(2) / 4 < rho_2 (x_i, x_j)$ ($2 <= i < j$), so $kappa_i^1 = 1$ for all $i >= 2$: $G_1$ is a star with $|E_1| = n - 1$, whereas $P P^top$ induces a clique on ${2, dots, n}$, giving $|E_A| = binom(n, 2) = Theta(n^2)$.

  *(2)* In the functional digraph $G_"dir" = (V, {(i, kappa_i^1)})$, every vertex has out-degree $1$, so each weakly connected component $C$ contains exactly one directed cycle $i_1 -> dots -> i_m -> i_1$ ($m >= 2$) with directed in-trees entering it. Writing $d_t = rho(x_(i_t), x_(i_(t + 1)))$ (mod $m$), optimality of $kappa_(i_t)^1$ and symmetry yield $d_t <= rho(x_(i_t), x_(i_(t - 1))) = d_(t - 1)$. Chaining gives $d_1 = dots = d_m$; distinct pairwise distances force $m = 2$, a unique MNN pair ${u, v}$. Outside ${u, v}$, any path step $i -> j -> k$ ($k != i$) obeys $rho(x_j, x_k) < rho(x_j, x_i)$ strictly. Merging the twin arcs $(u, v), (v, u)$ into one undirected edge leaves $|C| - 1$ edges on $C$, so $C$ is a tree.

  *(3)* Since $(P^2)_(i i) = sum_j P_(i j) P_(j i) = bb(I)[kappa_(kappa_i^1)^1 = i]$, each MNN pair contributes $2$ to $"Tr"(P^2)$. By (2), components biject with MNN pairs, so $n_1 = 1/2 "Tr"(P^2)$. Because $G_1$ is an $n_1$-tree forest on $n$ vertices, $|E_1| = n - n_1$.

  *(4)* For each $i in V$, ${i, kappa_i^1}$ is the unique lightest edge crossing the singleton cut $({i}, V without {i})$; by the MST cut property, $E_1 subset.eq E_"MST"$, matching Boruvka's first round identically.
]

Consequently, union-find over the $n$ arcs $(i, kappa_i^1)$ computes the exact FINCH partition in $O(n alpha(n))$ time and $O(n)$ space without ever forming $P P^top$.

#subsec("B", "Deterministic Logarithmic Depth and Poisson Compression")

#theorem(title: "Deterministic Logarithmic Depth and Tightness", origin: "Original Proof")[
  For any finite dataset $X$ ($n >= 2$, even with ties provided $kappa_a^1 != a$), FINCH's cluster counts satisfy
  $ n_ell <= floor(n_(ell - 1) / 2) <= floor(n / 2^ell) quad (ell >= 1), $
  terminating strictly at $n_L = 1$ in $L <= floor(log_2 n)$ steps. This bound is tight: for every $m >= 1$ and $n = 2^m$, there exists $X subset.eq bb(R)$ attaining $n_ell = 2^(m - ell)$ and $L = log_2 n$.
]

#proof[
  Whenever $n_(ell - 1) >= 2$, every vertex $a$ has an incident edge ${a, kappa_a^1}$ ($kappa_a^1 != a$) in $G_1^((ell))$, so every component satisfies $|C_(r, ell)| >= 2$. Hence $n_(ell - 1) = sum_(r = 1)^(n_ell) |C_(r, ell)| >= 2 n_ell$, giving $n_ell <= floor(n_(ell - 1) / 2) < n_(ell - 1)$. Using $floor(floor(x / a) / b) = floor(x / (a b))$, induction yields $n_ell <= floor(n / 2^ell)$. At $ell^* = floor(log_2 n)$, $1 <= n / 2^(ell^*) < 2$, forcing $n_(ell^*) = 1$ and $L <= floor(log_2 n)$.

  For tightness, map $b = (b_1, dots, b_m) in {0, 1}^m$ to the base-4 Cantor set
  $ x(b) = sum_(t = 1)^m b_t dot 4^(t - 1). $
  Pairs differing only in $b_1$ have distance $1$; pairs differing at maximal index $r >= 2$ satisfy
  $ |x(b) - x(b')| >= 4^(r - 1) - sum_(t = 1)^(r - 1) 4^(t - 1) = (2 dot 4^(r - 1) + 1) / 3 >= 3 > 1. $
  Thus $G_1^((1))$ is a perfect matching of $n_1 = 2^(m - 1)$ MNN pairs whose centroids
  $ mu^((1)) (b_2, dots, b_m) = 1 / 2 + 4 sum_(t = 1)^(m - 1) b_(t + 1) dot 4^(t - 1) $
  form an exact positive affine image of the $(m - 1)$-bit set. Induction yields $n_ell = 2^(m - ell)$ for all $ell <= m$ and $L = log_2 n$.
]

#proposition(title: "Expected Poisson Compression Ratio", origin: "Original Derivation")[
  On a homogeneous Poisson point process $Phi$ of intensity $lambda > 0$ in $(bb(R)^d, rho_2)$, Theorem 1 implies that FINCH's expected level-1 cluster-to-point ratio equals
  $ gamma_d = 1 / (2 (2 - I_d)), quad I_d = ("Vol"(B_d (0, 1) inter B_d (e_1, 1))) / ("Vol"(B_d (0, 1))). $
  Specifically, $gamma_1 = 1/3$, $gamma_2 = (8/3 + sqrt(3) / pi)^(-1) approx 0.3108$, and for all $d >= 1$,
  $ 1 / 4 < gamma_d <= 1 / (4 - 2 (sqrt(3) / 2)^d) -> 1 / 4 quad (d -> infinity). $
]

#proof[
  By Theorem 1, $gamma_d = p_d / 2$ where $p_d = "Pr"(kappa_(kappa_0^1)^1 = 0)$ for a typical point at $0$. By Slivnyak--Mecke, $Phi without {0}$ is an independent Poisson process of intensity $lambda$. Let $R = norm(Y)_2$ be the distance to $0$'s nearest neighbor $Y$, and $omega_d = "Vol"(B_d (0, 1))$; then $U = lambda omega_d R^d tilde.op "Exp"(1)$. Conditioned on $Y = y$ ($norm(y)_2 = r$), $B_d (0, r)$ is empty, so $kappa_y^1 = 0$ iff $B_d (y, r) without B_d (0, r)$---of volume $(1 - I_d) omega_d r^d$---is empty. Integrating over $U$ yields
  $ p_d = integral_0^infinity e^(-(1 - I_d) u) e^(-u) dif u = 1 / (2 - I_d), $
  giving $gamma_d = 1 / (2(2 - I_d))$. Direct lens areas give $I_1 = 1/2 => gamma_1 = 1/3$ and $I_2 = 2/3 - sqrt(3) / (2 pi) => gamma_2 approx 0.3108$. For any $z in B_d (0, 1) inter B_d (e_1, 1)$, the parallelogram identity gives
  $ norm(z - 1 / 2 e_1)_2^2 &= 1 / 2 norm(z)_2^2 + 1 / 2 norm(z - e_1)_2^2 - 1 / 4 norm(e_1)_2^2 \
                             &<= 1 / 2 + 1 / 2 - 1 / 4 = 3 / 4. $
  Hence $B_d (0, 1) inter B_d (e_1, 1) subset.eq B_d (1/2 e_1, sqrt(3)/2)$, implying $0 < I_d <= (sqrt(3)/2)^d$ and yielding the exponential convergence to $1/4$.
]

Proposition 1 explains why FINCH compresses typical moderate-to-high-$d$ datasets by $approx 4 times$--$6 times$ per round, yielding only $approx log_4 n$ partitions even for millions of points @sarfraz2019finch.

#subsec("C", "Invariances, Pure Recovery, and Asynchronous Collapse")

#theorem(title: "Invariances, Pure-Prefix Recovery, and Asynchronous Collapse", origin: "Original Counterexample")[
  Let $X = {x_1, dots, x_n} subset.eq bb(R)^d$ have distinct pairwise distances.
  + *Monotone Level-1 Invariance:* For any strictly increasing $g : [0, infinity) -> [0, infinity)$, $cal(P)_1 (X, g compose rho) = cal(P)_1 (X, rho)$.
  + *Similarity & Order Invariance:* Under $rho_2$, ${cal(P)_1, dots, cal(P)_L}$ is invariant to all similarities $T(x) = s Q x + t$ ($s > 0$, $Q in O(d)$, $t in bb(R)^d$) and input permutations.
  + *Pure-Prefix Recovery:* Let $X = union.big_(c = 1)^K X_c$ ($|X_c| >= 2$) satisfy
    $ Delta = max_c "diam"(X_c) < delta = min_(c != c') "dist"("conv"(X_c), "conv"(X_(c'))). $
    Under $rho_2$, whenever every $X_c$ retains $m_c^((ell - 1)) >= 2$ pure subclusters, $cal(P)_ell$ is strictly pure; if all classes share collapse depth $ell^*$, then $cal(P)_(ell^*) = {X_1, dots, X_K}$.
  + *Asynchronous Collapse Counterexample:* For every $delta / Delta > 1$, there exists a 6-point 2-cluster set $X_1 union X_2 subset.eq bb(R)$ with $"diam"(X_c) = Delta$ and $"dist"(X_1, X_2) = delta$ for which no raw FINCH partition equals ${X_1, X_2}$.
]

#proof[
  *(1)* Strictly increasing $g$ preserves all pairwise distance orderings and thus every $kappa_i^1$ and $G_1$.

  *(2)* Linearity of weighted means gives $mu_(r, ell) (T(X)) = s Q mu_(r, ell) (X) + t$, and orthogonality yields $rho_2 (T(u), T(v)) = s rho_2 (u, v)$ at every level $ell$, preserving all 1-NN rankings inductively; distinct distances remove order sensitivity.

  *(3)* Any pure subcluster $S subset.eq X_c$ has centroid in $"conv"(X_c)$. When $m_c^((ell - 1)) >= 2$, every subcluster in $X_c$ has a sibling within distance $<= "diam"("conv"(X_c)) <= Delta$, whereas any foreign centroid lies at distance $>= delta > Delta$. Thus all level-$ell$ links remain intra-class, giving exact recovery if all classes finish contracting at the same $ell^*$.

  *(4)* For any $M > 2$ ($Delta = 1$, $delta = M - 1$; rescale for arbitrary $Delta > 0$), set
  $ X_1 = {0, 1}, quad X_2 = {M, M + 0.1, M + 0.9, M + 1}. $
  At level 1, pairs $0 <-> 1$, $M <-> M + 0.1$, and $M + 0.9 <-> M + 1$ yield $n_1 = 3$ centroids $mu_1 = 0.5$, $mu_2 = M + 0.05$, and $mu_3 = M + 0.95$: $X_1$ has collapsed ($m_1^((1)) = 1$) while $X_2$ remains split ($m_2^((1)) = 2$). At level 2, $mu_2 <-> mu_3$ (dist $0.9$), while $mu_1$ has no intra-class neighbor left and links to $mu_2$ (dist $M - 0.45$), forming one component $mu_1 -> mu_2 <-> mu_3$ ($n_2 = 1$). Hence FINCH outputs cardinalities $(6, 3, 1)$, bypassing ${X_1, X_2}$ for all $M > 2$.
]

Theorem 3(4) explains why FINCH requires optional target-$K$ post-merging @sarfraz2019finch whenever true clusters differ in internal hierarchical depth, echoing Kleinberg's scale-consistency tension @kleinberg2002impossibility.

#subsec("D", "Amortized Complexity and High-Dimensional Concentration")

#theorem(title: "Amortized Complexity and Contrast Decay", origin: "Classical — Cited")[
  Let $T_"1NN" (m, d) <= c_d m^beta (log_2 m)^gamma$ ($beta >= 1, gamma >= 0$) bound 1-NN computation for $m$ points in $bb(R)^d$ @sarfraz2019finch. Then:
  + With union-find on $G_1^((ell))$ and weighted centroid recursion, FINCH's total runtime across all $L <= floor(log_2 n)$ levels satisfies
    $ T_"total" (n, d) <= (2^beta) / (2^beta - 1) T_"1NN" (n, d) + O(n d + n alpha(n)). $
    Exact brute-force search ($beta = 2, gamma = 0$) incurs $<= 4/3$ overhead over one 1-NN pass ($O(n^2 d)$), while tree / ANN indexing ($beta = 1, gamma = 1$) runs in $O(n log n)$ time @sarfraz2019finch.
  + Under Beyer et al.'s contrast-loss regime @beyer1999nearest, where $"Var"(rho(X_i, X_j)) / (bb(E)[rho(X_i, X_j)])^2 -> 0$ as $d -> infinity$, $D_"max"^((d)) / D_"min"^((d)) -> 1$ in probability, degrading spatial trees toward $Theta(n^2 d)$ scans and producing fragile hub-dominated 1-NN graphs.
]

#proof[
  *(1)* At level $ell in {0, dots, L - 1}$, FINCH runs 1-NN in $T_"1NN" (n_ell, d)$, union-find in $O(n_ell alpha(n_ell))$, and weighted centroid updates in $O(n_ell d)$. By Theorem 2, $n_ell <= n / 2^ell$ and $log_2 n_ell <= log_2 n$, so
  $ sum_(ell = 0)^(L - 1) T_"1NN" (n_ell, d) &<= c_d n^beta (log_2 n)^gamma sum_(ell = 0)^infinity 2^(-beta ell) \
                                             &= (2^beta) / (2^beta - 1) c_d n^beta (log_2 n)^gamma, $
  while $sum_(ell = 0)^(L - 1) n_ell <= 2 n$ bounds total centroid and union-find work by $O(n d + n alpha(n))$.

  *(2)* By Chebyshev's inequality @beyer1999nearest, vanishing relative variance concentrates all pairwise distances near their mean so $D_"max"^((d)) / D_"min"^((d)) -> 1$ in probability. Tree pruning fails (query balls overlap almost all partitions) and tiny coordinate noise flips 1-NN links $kappa_i^1$, creating spurious star hubs that merge unrelated classes.
]

#sec("III", "Complexity, Geometry and Operational Profile")

FINCH's runtime is dominated by the level-1 nearest-neighbor pass (Theorem 4). For $n <= 5 times 10^4$, tiled BLAS-3 GEMM blocks compute exact cosine or Euclidean 1-NNs in $O(n^2 d)$ flops using only $O(n d + B n)$ memory (tile size $B$), avoiding $O(n^2)$ distance storage. At million scale, KD-trees (low $d$) or ANN graphs such as PyNNDescent / HNSW (high $d$) yield $O(n log n)$ scaling @sarfraz2019finch. Replacing explicit $P P^top$ multiplication with union-find on $G_1$ (Theorem 1) enforces strict $O(n alpha(n))$ time and $O(n)$ memory regardless of hubness.

Geometrically, FINCH pairs *two distinct regimes*: level 1 is purely local and ordinal (Theorem 3(1)), following curved manifolds and varying densities, whereas levels $ell >= 2$ contract arithmetic centroids, shifting inductive bias toward convex centroid-linkage regions (UPGMC) @lance1967general. Lacking a noise label, every outlier is absorbed at level 1, causing a $1/n$ centroid breakdown point.

#figure(
  table(
    columns: (1.1fr, 1.9fr),
    align: (left, left),
    [Property], [Characterization],
    [Time Complexity], [$O(n log n)$ ANN / tree @sarfraz2019finch; $O(n^2 d)$ exact ($<= 4/3$ total-hierarchy factor, Thm. 4)],
    [Space Complexity], [$O(n d)$ data + centroids; $O(n)$ union-find on $G_1$ (avoids $Theta(n^2)$ $P P^top$, Thm. 1)],
    [Hierarchy Depth], [$L <= floor(log_2 n)$ worst case (Thm. 2); $approx log_4 n$ typical ($gamma_d -> 1/4$, Prop. 1)],
    [Hyperparameters], [Zero required (metric $rho$ only; optional target-$K$ post-merge @sarfraz2019finch)],
    [Invariances], [Level-1 monotone metric; full-hierarchy similarity & permutation invariance (Thm. 3)],
    [Cluster Geometry], [Hybrid: local manifold-following at level 1; convex centroid-based at levels $ell >= 2$],
    [Outlier Robustness], [Low: $1/n$ centroid breakdown; no noise class (outliers merged at level 1)],
    [Determinism], [Exact under deterministic 1-NN; stochastic under randomized ANN indexes],
    [Parallelism], [Embarrassingly parallel 1-NN (SIMD / GPU GEMM); linear-time union-find],
  ),
  caption: [Operational and theoretical profile of FINCH (Alg. 80).],
)

#sec("IV", "Engineering and Scientific Applications")

#subsec("A", "Primary Application Domains and Mathematical Fit")

+ *Large-Scale Visual and Face Identity Clustering:* Deep metric encoders embed images on $bb(S)^(d - 1)$ where same-identity samples form tight angular neighborhoods with unknown, imbalanced class counts $K$ @sarfraz2019finch, @yang2016jule, @caron2018deepcluster. Cosine FINCH links local pose variations at levels 1--2 without tuning $K$ or density thresholds, scaling to millions of embeddings via ANN @sarfraz2019finch, @vangansbeke2020scan.
+ *Video Action Segmentation and Shot Discovery:* Frame-level descriptors from untrimmed video streams exhibit highly variable action counts across recordings @sarfraz2019finch. FINCH's $O(n log n)$ speed and $3$--$6$ level hierarchy jointly capture fine sub-actions (levels 1--2) and coarse scenes (levels 3--4) with zero per-video tuning.
+ *Single-Cell Transcriptomics Metacell Coarsening:* Large scRNA-seq atlases ($10^5$--$10^6$ cells) benefit from fast, parameter-free compression into homogeneous metacells. Level-1 FINCH reduces cardinality by $approx 4 times$ (Proposition 1) via purely local 1-NN links, preserving ultra-rare states of $2$--$5$ cells often discarded by density cutoffs @ester1996dbscan or fixed-$k$ methods @arthur2007kmeanspp.
+ *Speaker Diarization and Acoustic Indexing:* Normalized speaker embeddings (x-vectors / d-vectors) have unknown speaker cardinality and strict latency budgets. Cosine FINCH's similarity invariance (Theorem 3(2)) and sub-millisecond runtime partition recordings across speakers without BIC model-selection sweeps @schwarz1978bic.
+ *Parameter-Free Coreset Summarization:* Analogous to BIRCH micro-clusters @zhang1996birch, FINCH's level-1 / level-2 weighted centroids $(mu_(r, ell), w_(r, ell))$ preserve exact Euclidean variance statistics, cutting sample size by $approx 4^ell$ before expensive $O(n^2)$--$O(n^3)$ spectral @ng2002spectral or Ward @ward1963hierarchical refinement.

#subsec("B", "Failure Modes and Domain Suitability Matrix")

Four structural failure modes limit FINCH: (i) *asynchronous multi-scale collapse* (Theorem 3(4)), skipping ground-truth partitions when clusters have unequal internal depths; (ii) *outlier absorption*, corrupting level-$>= 2$ centroids ($1/n$ breakdown); (iii) *high-dimensional contrast loss and hubness* (Theorem 4(2), @beyer1999nearest), collapsing classes via star hubs; and (iv) *macro-scale non-convexity*, where concentric rings or spirals share overlapping centroids at levels $ell >= 2$.

#figure(
  table(
    columns: (1.0fr, 1.2fr, 0.8fr),
    align: (left, left, left),
    [Domain], [Mathematical Fit], [Caution / Alternative],
    [Deep Vision & Faces], [High: cosine 1-NN tracks identity manifolds; unknown $K$], [Use ANN + $K$-merge; compare SCAN @vangansbeke2020scan],
    [Video & Speech Streams], [High: sub-ms multi-scale hierarchy; zero tuning], [Post-merge for exact $K$ under varying scales],
    [scRNA-seq Metacells], [High at levels 1--2: $approx 4 times$ local pure compression], [Run PCA first; use Leiden @traag2019leiden globally],
    [Noisy Spatial / LiDAR], [Poor: $1/n$ outlier absorption skews centroids], [Prefer HDBSCAN @campello2013hdbscan or RNN-DBSCAN @bryant2018rnndbscan],
    [Concentric / Spirals], [Mixed: level 1 chains curves; levels $>= 2$ average centroids], [Prefer Spectral @ng2002spectral or SLINK @sibson1973slink],
    [Raw High-$d$ Tabular], [Low: contrast decay & hubness @beyer1999nearest], [Reduce dimension first or use SNN @ertoz2003snn],
  ),
  caption: [Domain suitability matrix and failure-mode alternatives for FINCH.],
)

#sec("V", "Structural Relations and Empirical Behavior")

Within Cat. III, FINCH bridges single-linkage MST clustering (SLINK @sibson1973slink, @gower1969minimum) and centroid linkage (UPGMC @lance1967general, @mullner2011modern), while sharing reciprocal-neighbor principles with SNN @ertoz2003snn, CHAMELEON @karypis1999chameleon, RNN-DBSCAN @bryant2018rnndbscan, PERCH @kobren2017perch, and SCC @monath2021scc. Our final proposition delineates exactly when FINCH level 1 diverges from SLINK and proves its $1/n$ centroid breakdown point.

#proposition(title: "Local Adaptivity vs. Single-Linkage and Outlier Breakdown", origin: "Original Proof")[
  Let $(X, rho)$ have distinct distances, unique MST edge set $E_"MST"$, and level-1 FINCH forest $E_1 subset.eq E_"MST"$ (Theorem 1). Then:
  + $cal(P)_1$ equals a global Single-Linkage (SLINK @sibson1973slink) dendrogram cut if and only if
    $ max_(e in E_1) w(e) < min_(e' in E_"MST" without E_1) w(e'). $
  + In general, $cal(P)_1$ need not appear anywhere in the SLINK dendrogram, and level-$>= 2$ centroids have breakdown point $1/n$.
]

#proof[
  *(1)* Every SLINK cut at threshold $tau$ retains precisely ${e in E_"MST" : w(e) <= tau}$ @gower1969minimum. Since $E_1$ and any cut are both subforests of $E_"MST"$, their partitions coincide iff their edge sets coincide, which holds for some $tau$ iff every edge in $E_1$ is strictly shorter than every edge in $E_"MST" without E_1$.

  *(2)* Take $X = {0, 1, 3, 3.5, 10} subset.eq bb(R)$ ($n = 5$). Sorted MST edges are ${3, 3.5}$ ($0.5$), ${0, 1}$ ($1.0$), ${1, 3}$ ($2.0$), and ${3.5, 10}$ ($6.5$), so SLINK only produces ${{0}, {1}, {3, 3.5}, {10}}$, ${{0, 1}, {3, 3.5}, {10}}$, and ${{0, 1, 3, 3.5}, {10}}$. Because $kappa^1 (10) = 3.5$, FINCH yields $cal(P)_1 = {{0, 1}, {3, 3.5, 10}}$, which never appears in SLINK ($6.5 > 2.0$). Sending $10 -> M -> infinity$ preserves $kappa^1 (M) = 3.5$ and drives the centroid $(6.5 + M) / 3 -> infinity$, proving breakdown point $1/n$.
]

*Empirical Benchmark Analysis.* On the ClusterOverview 4-dataset benchmark evaluated by Silhouette @rousseeuw1987silhouettes and Adjusted Rand Index (ARI) @hubert1985comparing, FINCH closely confirms Sections II--V:
- *Iris* (`Sil = 0.586`, `ARI = 0.567`): _Setosa_ isolates cleanly at level 1, while overlapping _Versicolor_ / _Virginica_ points trigger cross-class 1-NN links and asynchronous collapse (Theorem 3(4)), producing strong Silhouette (`0.586`) with moderate ARI (`0.567`).
- *Customers* (`Sil = 0.359`): Skewed high-spend outliers are absorbed into adjacent clusters at level 1 without noise filtering (Proposition 2), moderating geometric cohesion (`0.359`).
- *Multi-Shape* (`Sil = 0.712`, `ARI = N/A`): Locally adaptive 1-NN chaining paired with compact centroid merging achieves FINCH's highest Silhouette (`0.712`).
- *High-D* (`Sil = 0.185`, `ARI = N/A`): Contrast decay and hubness (Theorem 4(2), @beyer1999nearest) plus rapid $gamma_d -> 1/4$ compression (Proposition 1) depress Silhouette to `0.185`.
- *Runtime* (`Mean = 0.569 ms`): Confirming Theorems 2 and 4, parallel logarithmic contraction ($L <= floor(log_2 n)$) delivers sub-millisecond execution (`0.569 ms`), ranking among the fastest parameter-free hierarchical algorithms in the catalog.
