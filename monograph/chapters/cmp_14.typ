#import "../template/ieee.typ": *

#cmp_entry(
  roman: "XIV",
  title: "Practitioner's Engineering Applications Handbook and Algorithm Selection Architecture",
  scope_label: "Domain-to-Algorithm Mathematical Selection Framework Across 10 Scientific & Industrial Verticals",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

Comparison Chapters I--XIII have established the complete theoretical and empirical foundation of this 200-algorithm compendium: six unifying variational, spectral, and topological lenses (Chapter I), fifteen exact cross-family reductions (Chapter II), intra-category trade-offs across Categories I--VIII (Chapters III--X), axiomatic impossibility and concentration meta-theorems (Chapter XI), the 200-row ten-dimension reference compendium (Chapter XII), and controlled multi-geometry benchmarks in the ClusterOverview engine (Chapter XIII). Yet an urgent engineering question governs every real deployment: _given only observable data invariants—cardinality, dimension, metric geometry, sparsity, contamination rate, streaming memory bounds, and auditability mandates—how does a practitioner map an unfamiliar workload onto a provably sound subset of Algorithms 1--200 without blind trial-and-error?_

By Kleinberg's impossibility theorem @kleinberg2002impossibility and unsupervised No-Free-Lunch geometry @jain2010data, @vonluxburg2007tutorial, no single algorithm dominates all data-generating processes. Algorithm selection is nevertheless governed by exact *dual elimination and constructive matching*: first pruning algorithm families whose mathematical preconditions conflict with observed data invariants, then assembling a compact cross-paradigm ensemble whose inductive biases span the surviving structural hypotheses @xu2005survey, @aggarwal2013data, @strehl2002cluster.

This closing handbook delivers that engineering architecture in four parts. Section I formalizes the *Algorithm Selection Architecture*, proving a six-rule Master Elimination and Soundness Theorem (Theorem 1) and a four-stage selection funnel (Proposition 1 and Algorithm 1). Section II synthesizes domain pipelines, failure modes, and formal recovery laws across *ten major scientific and industrial verticals* (Propositions 2--7 and Theorems 2--5). Section III consolidates the complete mapping into two full-width Master Diagnostic Decision Matrices (Tables I--II). Section IV establishes a *Production Engineering Checklist* covering scaling invariance (Proposition 8), seven recurring validation pitfalls, and a formal cross-paradigm verification certificate (Theorem 6 and Table III) operationalized through the pure-Rust ClusterOverview engine.

// =============================================================================
// SECTION I: FORMAL ALGORITHM SELECTION ARCHITECTURE
// =============================================================================
#sec("I", "Formal Algorithm Selection Architecture and Decision Theorems")

To replace informal heuristics with verifiable engineering rules, we summarize a clustering workload by its observable and estimable *data-invariant signature* and characterize its admissible cone inside the 200-algorithm catalog $cal(A)_(200) := {1, dots, 200}$.

#definition(title: "Data-Invariant Signature and Admissible Algorithm Cone")[
  Let $cal(X) = {x_1, dots, x_n}$ be drawn from an unknown process $cal(P)$ on ambient domain $Omega$ ($bb(R)^d$, mixed product space $bb(R)^(d_r) times product_(q=1)^(d_c) cal(V)_q$, unit sphere $bb(S)^(d-1)$, graph $G = (V, E, W)$, sequence space, or two-mode matrix $bb(R)^(n times p)$). Its *workload signature* is the ten-tuple
  $ cal(I)(cal(X)) := ( & n, d, d_"int", cal(M), rho_"sp", \
    & eta_"out", cal(G), B_"mem", tau_"lat", iota_"int" ), $
  encoding: cardinality $n$ ($n = infinity$ for streams); ambient and intrinsic dimensions $d_"int" <= d$; native dissimilarity structure $cal(M)$; sparsity fraction $rho_"sp" in [0, 1]$; Huber outlier fraction $eta_"out" in [0, 1 / 2)$; cluster geometry $cal(G)$ (spherical, ellipsoidal, non-convex manifold, multi-density, overlapping, or union of subspaces); RAM and latency budgets $(B_"mem", tau_"lat")$; and interpretability constraint $iota_"int" in {"strict", "moderate", "blackbox"}$. An algorithm $A in cal(A)_(200)$ lies in the *$epsilon$-admissible cone* $cal(C)_epsilon (cal(I))$ if $A$ meets $(B_"mem", tau_"lat", iota_"int")$ and attains consistent or $(1 + epsilon)$-bounded risk on the model class of $(cal(M), d_"int", rho_"sp", eta_"out", cal(G))$.
]

#subsec("A", "Master Algorithm Selection Soundness and Elimination Theorem")

Rigorous selection starts with *negative pruning*: eliminating families in $cal(A)_(200)$ whose mathematical preconditions are provably violated by $cal(I)(cal(X))$.

#theorem(title: "Master Algorithm Selection Soundness and Elimination Theorem", origin: "Original Unified Proof")[
  For any workload signature $cal(I)(cal(X))$ with $k >= 2$ true clusters, six sharp elimination and soundness laws govern Categories I--VIII (Algorithms 1--200):
  1. *(Memory & Streaming Rule):* With scalar width $w_"flt"$ ($8$ bytes for `f64`), if $n (n - 1) w_"flt" / 2 > B_"mem"$, every Tier 4--5 method materializing an $O(n^2)$ pairwise matrix (PAM \#4, dense Spectral \#103--104, exact linkages \#56--62, \#73, \#75, CSPA \#171) is memory-infeasible. For $B_"mem" = o(n d)$, single-pass $(1 + epsilon)$-approximation to $k$-Means is achieved by merge-and-reduce coresets (StreamKM++ \#195, BICO \#78) in $O(k d epsilon^(-2) log(k / delta) log^2 n)$ space @harpeled2004coresets, @ackermann2012streamkmpp, @fichtenberger2013bico.
  2. *(High-Dimensional Concentration Rule):* Let $cal(X) subset.eq bb(R)^d$ combine $d_"int"$ signal coordinates with $d - d_"int"$ independent noise dimensions of variance $sigma_"bg"^2 > 0$. If separation $Delta_"sig" := min_(a != b) norm(mu_a - mu_b)_2$ satisfies
     $ Delta_"sig"^2 = o(sigma_"bg"^2 sqrt((d - d_"int") log n)), $
     ambient $L_2$ contrast vanishes @beyer1999nearest, causing isotropic Voronoi partitioners (Cat. I, Algs. 1--11, 13--25), ambient fixed-$epsilon$ density extractors (Algs. 31, 33--35, 38, 47--50), and spatial trees ($d_"int" > 15$) to collapse. Admissible survivors require directional $bb(S)^(d-1)$ geometry (Algs. 12, 89), SNN rank overlap (Algs. 36, 52), subspace projection (Cat. VI, Algs. 134--150), sparse topic/matrix models (Algs. 83, 91--92, 97, 100--101), or deep latent encoders (Cat. VIII, Algs. 181--192).
  3. *(Contamination Breakdown Rule):* Under Huber contamination $P_eta = (1 - eta) P_0 + eta Q$ ($eta > 0$), arithmetic-mean and Gaussian-MLE estimators ($k$-Means \#1--3, \#9--14, \#17--23; Ward \#59; GMM \#81, \#93; FCM \#151) have breakdown point $epsilon^* = 1 / n$ and unbounded drift $sup_Q norm(hat(mu)_j - mu_j)_2 = infinity$. Bounded error strictly requires $alpha$-trimming with $alpha >= eta$ (Trimmed K-Means \#26 @cuesta1997trimmed), $L_1$/medoid estimators (PAM \#4, FasterPAM \#24, K-Medians \#15), density-excision noise sets (Algs. 31--33, 53, 194), or redescending Student-$t$ mixtures (TMM \#88 @peel2000robust).
  4. *(Non-Convexity & Density-Ratio Rule):* If $"conv"(C_a) inter C_b != emptyset$ for some $a != b$, *every* Bregman centroid partitioner (Cat. I, isotropic GMM \#81) incurs strictly positive Bayes error. Further, if two clusters have peak densities $f_1 > f_"val" > f_2 > 0$ separated by saddle $f_"val"$, *no* single global level-set threshold (hence no fixed-$(epsilon, "minPts")$ DBSCAN \#31) can separate the dense pair while preserving the sparse cluster—mandating multi-scale density/reachability (HDBSCAN \#32, OPTICS \#33, VDBSCAN \#50, DeBaCl \#46, SNN \#36) or normalized spectral/graph cuts (Algs. 103--104, 107--108, 113).
  5. *(Domain-Metric Invariance Rule):* Embedding non-Euclidean spaces (unit spheres, count simplices, mixed numeric-categorical tables, warped time series, or sparse graphs) into unweighted $L_2$ injects nuisance distortion; invariance restricts candidates to native operators: cosine/vMF (Algs. 12, 89), multinomial/Dirichlet/Poisson mixtures (Algs. 83, 91, 95--97, 100--101), mixed medoids/prototypes (Algs. 4, 7--8, 24, 67, 90), elastic/shape barycenters (Algs. 162--164), or graph/SBM models (Algs. 103--132).
  6. *(Interpretability Rule):* When auditability demands $iota_"int" = "strict"$ (explicit coordinate attribution), nonlinear deep encoders (Algs. 181--192) and implicit kernels (Algs. 16, 155) are disqualified; exact coordinate provenance restricts $cal(C)_epsilon (cal(I))$ to axis-aligned subspaces (Algs. 134, 136--137, 145, 148--149), topic/mixture models (Algs. 81, 83, 88--94, 97, 100), medoids (Algs. 4, 8, 24, 106), biclusters (Algs. 157--161), and dendrograms (Algs. 56--63, 73, 75).
]

#proof[
  *(Rule 1)* Storing the upper triangle of $D in bb(R)^(n times n)$ requires $(n (n - 1) / 2) w_"flt"$ bytes; exceeding $B_"mem"$ precludes exact $O(n^2)$ materialization. Conversely, Bentley--Saxe merge-and-reduce in StreamKM++ (\#195 @ackermann2012streamkmpp) maintains $L = ceil(log_2 (n / m))$ buffers of size $m = O(k d epsilon^(-2) log(k / delta) log n)$ with per-level distortion $1 + epsilon / (2 L)$, bounding total distortion across $L$ merges by
  $ (1 + epsilon / (2 L))^L <= exp(epsilon / 2) <= 1 + epsilon quad (0 < epsilon <= 1) $
  in strictly $O(L m d) = o(n d)$ memory.

  *(Rule 2)* Splitting $x_i = (x_i^((S)), x_i^((N))) in bb(R)^(d_"int") times bb(R)^(d_"noise")$ ($d_"noise" = d - d_"int"$),
  $ norm(x_i - x_j)_2^2 &= norm(x_i^((S)) - x_j^((S)))_2^2 \
    &quad + norm(x_i^((N)) - x_j^((N)))_2^2. $
  The noise term $Z_(i j) := norm(x_i^((N)) - x_j^((N)))_2^2$ has mean $2 d_"noise" sigma_"bg"^2$ and extremal fluctuations $Theta(sigma_"bg"^2 sqrt(d_"noise" log n))$ across pairs. When $Delta_"sig"^2 = o(sigma_"bg"^2 sqrt(d_"noise" log n))$, random noise swamps inter-cluster signal, driving nearest-neighbor purity to chance $1 / k$. Projecting onto $S$ (Cat. VI) removes $Z_(i j)$ completely; $bb(S)^(d-1)$ cosine models cancel norm variance; and SNN rank overlap retains positive correlation under structured noise.

  *(Rule 3)* Appending one adversarial point $x_"adv" = R bold(e)_1$ ($R -> infinity$) either pulls an unweighted centroid by $norm(hat(mu)_j)_2 >= R / n -> infinity$ or forces a singleton cluster that collapses two genuine clusters (and sends GMM variance $sigma_j^2 -> 0$ with unbounded likelihood), proving $epsilon^* = 1 / n$. In Trimmed K-Means (\#26 @cuesta1997trimmed) with $alpha >= eta$, the $ceil(alpha n)$ largest residuals—including all $x_"adv"$ outside clean radius—are dropped before mean updates; in PAM/K-Medians (\#4, \#15, \#24), unit-bounded $L_1$ subgradients prevent drift up to half the smallest cluster; in TMM (\#88 @peel2000robust), $u_(i j) norm(x_i - mu_j)_2 -> 0$ as $norm(x_i)_2 -> infinity$.

  *(Rule 4)* By Banerjee et al. @banerjee2005clustering, Bregman Voronoi cells are intersections of half-spaces, hence convex; if $x^* in "conv"(C_a) inter C_b$, any convex cell covering $C_a$ misclassifies $x^*$. For density ratios, separating two peaks across saddle $f_"val"$ requires upper level-set threshold $lambda > f_"val"$, whereas retaining a sparse peak $f_2 < f_"val"$ requires $lambda <= f_2$; since $(f_"val", f_2] = emptyset$, no global threshold succeeds. HDBSCAN (\#32 @campello2013hdbscan) escapes via branch-local excess-of-mass cuts at distinct heights $lambda_1 != lambda_2$.

  *(Rule 5)* Rescaling identical unit vectors $u = v in bb(S)^(d-1)$ by document lengths $L_u != L_v$ gives $norm(L_u u - L_v v)_2^2 = (L_u - L_v)^2 > 0$; shifting $sin(omega t)$ by $pi$ maximizes $L_2$ despite zero SBD (\#162) or Soft-DTW (\#164) dissimilarity. Metric-native operators quotient out nuisance transformations identically.

  *(Rule 6)* Pulling back a latent linear Voronoi cell $V_j$ through an $L$-layer nonlinear encoder $f_theta$ yields a non-convex piecewise-algebraic set $f_theta^(-1) (V_j)$ lacking coordinate-level thresholds, whereas axis-aligned subspaces, topic posteriors, medoids, and ultrametric trees provide explicit coordinate certificates.
]

#subsec("B", "Four-Stage Lexicographic Selection Funnel")

Theorem 1 compiles directly into a constant-depth four-stage engineering pipeline reducing $cal(A)_(200)$ to a diverse shortlist of $3$--$6$ algorithms.

#proposition(title: "Completeness and Complexity of the Four-Stage Selection Funnel", origin: "Original Proof")[
  Let $cal(S)_1 supset.eq cal(S)_2 supset.eq cal(S)_3 supset.eq cal(S)_4 subset.eq cal(A)_(200)$ be filtered sequentially by: (Stage 1) Metric Compatibility (Rule 5); (Stage 2) Scale & Streaming Memory Gate (Rule 1); (Stage 3) Dimension, Geometry, Noise & Audit Gate (Rules 2, 3, 4, 6); and (Stage 4) Cross-Paradigm Diversification selecting $1$--$2$ Pareto-efficient representatives per surviving category. Then: (i) every pruned method violates a necessary feasibility or consistency condition of Theorem 1; (ii) for every realizable $cal(I)(cal(X))$, $3 <= |cal(S)_4| <= 6$ spans $>= 2$ distinct paradigms (whenever $B_"mem" = Omega(n d)$); and (iii) executing $cal(S)_4$ plus all pairwise ARI cross-checks adds only $O(n + k^2)$ overhead over running a single admissible candidate.
]

#proof[
  *(i)* Immediate from Theorem 1 Rules 1--6. *(ii)* Chapter XII Table I confirms all eight categories cover every invariant regime (32 Tier-1 streaming methods, 102 high-$d$ methods, 103 noise-robust methods, 102 auto-$k$ methods), ensuring non-empty multi-category survival capped at $|cal(S)_4| <= 6$. *(iii)* Building the $k times k$ contingency table for pairwise ARI @hubert1985comparing costs $O(n)$ per pair; across $binom(|cal(S)_4|, 2) <= 15$ pairs, total verification overhead is $O(n + k^2)$, dominated by the $Omega(n d)$ coordinate-reading lower bound (Chapter XII, Theorem 1).
]

#algbox(
  title: "Four-Stage Algorithm Selection Architecture (Select200)",
  inputs: [Signature $cal(I) = (n, d, d_"int", cal(M), rho_"sp", eta_"out", cal(G), B_"mem", tau_"lat", iota_"int")$],
  outputs: [Cross-paradigm shortlist $cal(S)_4 subset.eq cal(A)_(200)$ ($3 <= |cal(S)_4| <= 6$) + policy $Pi_"val"$],
  [
    1. *Stage 1 (Metric Gate):* Filter $cal(A)_(200) -> cal(S)_1$ by $cal(M)$: _Graph:_ Cat. V (\#103--132); _Text/Sphere:_ \#12, \#83--84, \#89, \#91--92, \#95--97, \#100--101, \#111; _Mixed:_ \#4--8, \#24--25, \#67, \#85--87, \#90, \#94; _Time Series:_ \#162--165, \#169; _2-Mode:_ \#123, \#157--161; _Continuous/Multi-View:_ Euclidean, Density, Subspace, Multi-View & Deep pools.
    2. *Stage 2 (Scale & RAM Gate):* If stream or $n d w_"flt" > B_"mem"$, retain Tier-1 streaming/coreset engines (\#3, \#23, \#64, \#78, \#193--197). Else if $n (n - 1) w_"flt" / 2 > B_"mem"$, eliminate dense $O(n^2)$ matrix methods (Tier 4--5).
    3. *Stage 3 (Geometry, High-$d$, Noise & Audit Gate):* If $d > 50, d_"int" << d$: require Subspace, Cosine/Topic, SNN, or Deep Latent (Rule 2). If $eta_"out" > 0.02$: require Trimmed (\#26), Medoid (\#4, \#24), Density (\#31--33, \#53), or TMM (\#88) (Rule 3). If non-convex/multi-density: enforce HDBSCAN/OPTICS/Spectral/Leiden (Rule 4). If $iota_"int" = "strict"$: exclude deep encoders (\#181--192) and implicit kernels (\#16, \#155) (Rule 6).
    4. *Stage 4 (Diversification & Verification Policy):* Pick $1$--$2$ top methods per surviving category ($>= 2$ paradigms) for $cal(S)_4$; bind $Pi_"val"$ prohibiting Silhouette-only selection on non-convex data and requiring ClusterOverview A-vs-B ARI verification.
  ],
)

// =============================================================================
// SECTION II: TEN SCIENTIFIC AND INDUSTRIAL ENGINEERING VERTICALS
// =============================================================================
#sec("II", "Engineering Deployment Across Ten Scientific and Industrial Verticals")

We now instantiate the Section I architecture across ten flagship scientific and industrial verticals, diagnosing domain failure modes, proving vertical-specific recovery laws, and specifying exact Tier-1 pipelines across Algorithms 1--200.

#subsec("A", "Vertical 1: Genomics, Transcriptomics, and Systems Biology")

*Data Invariants and Failure Modes:* Molecular biology spans three regimes: (i) *single-cell RNA-seq (scRNA-seq)* count matrices $X in bb(Z)_(>= 0)^(n times p)$ ($n in [10^4, 10^7]$, $p approx 2 times 10^4$ genes, sparsity $rho_"sp" in [0.85, 0.98]$ from technical dropouts, negative-binomial overdispersion, and continuous lineage manifolds); (ii) *bulk expression and drug-response matrices* where biological pathways co-regulate subsets of genes *only across subsets of conditions* (overlapping two-mode biclusters); and (iii) *phylogenetic sequence-distance matrices* exhibiting lineage-specific mutation-rate heterogeneity. Raw $k$-Means (\#1) fails on scRNA-seq via high-$d$ concentration (Rule 2), library-size confounding (Rule 5), and Voronoi absorption of rare $< 0.5\%$ stem-cell populations.

*Canonical Algorithm Architecture:*
- *scRNA-seq Cell Atlases:* After library-size normalization, log1p/variance stabilization, highly variable gene selection ($p' approx 2{,}000$), and truncated PCA ($d_"pca" in [20, 50]$), pipelines construct a *Shared Nearest Neighbor (SNN)* graph (Alg. 36 @ertoz2003snn) partitioned via the *Leiden algorithm* (Alg. 108 @traag2019leiden), replacing *Louvain* (Alg. 107 @blondel2008louvain) to guarantee connected cell communities. Where ambient-RNA debris or transitional doublets abound, *HDBSCAN* (Alg. 32 @campello2013hdbscan) on diffusion/PCA embeddings cleanly excises noise and resolves rare subpopulations without pre-specifying $k$.
- *Expression Biclustering & Oncology Subtyping:* *Cheng–Church* (Alg. 157 @cheng2000biclustering) extracts additive coherent gene-condition modules minimizing Mean Squared Residue $"MSR"(I, J) <= delta$; *FABIA* (Alg. 161 @hochreiter2010fabia) fits multiplicative overlapping sparse factor models $X = sum_(l=1)^p lambda_l z_l^top + Upsilon$ with Laplacian priors, complemented by *Kluger Spectral Biclustering* (\#123 @kluger2003spectral), *ISA* (\#160 @bergmann2003isa), *PLAID* (\#158 @lazzeroni2002plaid), and binary *Bimax* (\#159 @prelic2006bimax). For TCGA-scale cancer subtyping ($n in [80, 1500]$), *Consensus Clustering* (Monti et al., Alg. 170 @monti2003consensus) aggregates $B in [200, 1000]$ subsampled runs into co-association CDF curves to stabilize $k$.
- *Phylogenetics:* *UPGMA* (Alg. 58 @sokal1958statistical) builds rooted ultrametric trees under strict molecular clocks, whereas *Neighbor-Joining (NJ)* (Alg. 75 @saitou1987neighbor) reconstructs exact unrooted additive trees despite unequal evolutionary rates.

#proposition(title: "Phylogenetic Additive Consistency and Leiden Connectedness in Genomics", origin: "Classical — Cited")[
  1. *(UPGMA vs. Neighbor-Joining):* For an additive tree metric $d_T$ on $n$ taxa with shortest internal branch $e_"min" > 0$, UPGMA (Alg. 58) recovers $T$ iff $d_T$ is ultrametric. By contrast, Neighbor-Joining (Alg. 75 @saitou1987neighbor), minimizing the Saitou–Nei criterion
     $ Q(i, j) := (n - 2) d(i, j) - sum_(r=1)^n d(i, r) - sum_(r=1)^n d(j, r), $
     recovers the exact unrooted topology of *every* additive tree and tolerates perturbations $norm(hat(d) - d_T)_infinity < e_"min" / 2$.
  2. *(Leiden Connectedness vs. Louvain Pathology):* On sparse scRNA-seq SNN graphs, Louvain (Alg. 107 @blondel2008louvain) can output communities whose induced subgraphs $G[C]$ are disconnected, whereas Leiden's refinement step (Alg. 108 @traag2019leiden) guarantees every returned community is internally connected and $gamma$-separated.
]

#proof[
  *(1)* On the 4-taxon additive tree $((A, B), (C, D))$ with internal edge $e = 1$ and terminal branches $a = d = 10, b = c = 1$, exact distances satisfy $d(A, B) = d(C, D) = 11$ but $d(B, C) = 3$. UPGMA erroneously merges $(B, C)$ first (long-branch attraction). With row sums $R_A = R_D = 44, R_B = R_C = 26$, NJ yields
  $ Q(A, B) &= 2(11) - 44 - 26 = -48, \
    Q(B, C) &= 2(3) - 26 - 26 = -46 > -48, $
  selecting the true cherry $(A, B)$ strictly. Generally, the four-point condition bounds $Q(i, j) >= min_"cherries" Q(u, v) + 2 e_"min"$ for every non-cherry pair, yielding Atteson's $e_"min" / 2$ $L_infinity$-stability radius.

  *(2)* In Louvain, a bridge cell $v$ linking two cliques $K_a, K_b subset.eq C$ can migrate to $C'$ for positive modularity gain, leaving $C without {v}$ disconnected prior to irreversible supernode aggregation. Leiden refines each community $C$ from singletons before aggregation, merging only subsets $S subset.eq C$ satisfying $E(S, C without S) >= gamma |S|(|C| - |S|) > 0$; induction guarantees internal connectivity of every aggregated community @traag2019leiden.
]

#subsec("B", "Vertical 2: Computer Vision, Medical Imaging, and 3D LiDAR Robotics")

*Data Invariants and Failure Modes:* Visual and robotic sensing spans three substrates: (i) *2D RGB and 3D MRI/CT volumes* ($N in [10^5, 10^8]$ pixels/voxels) with spatial autocorrelation, partial-volume boundaries, and Rician noise; (ii) *autonomous-driving 3D LiDAR point clouds* ($N in [6 times 10^4, 3 times 10^5]$ points/scan in $bb(R)^3$ at 10 Hz) featuring radial $r^(-2)$ density falloff and $< 20$ ms latency budgets; and (iii) *large-scale natural image corpora* ($N in [10^5, 10^7]$) where raw pixel $L_2$ ignores pose, illumination, and semantic invariance. Global intensity $k$-Means (\#1) fragments noisy MRI masks into salt-and-pepper artifacts, while pixel-level clustering collapses on semantic recognition.

*Canonical Algorithm Architecture:*
- *Superpixels, Watershed & Fast Mode Seeking:* *SLIC* (Alg. 166 @achanta2012slic) generates compact superpixels in $O(N)$ time by restricting 5D CIELAB+$(x, y)$ updates to local $2 S times 2 S$ windows. *Watershed* (Alg. 167 @vincent1991watersheds) segments touching cells or anatomical contours on gradient relief in $O(N log N)$ time using marker-controlled flooding. For nonparametric color-spatial mode seeking without fixing $k$, *Mean Shift++* (Alg. 55 @jang2021meanshiftpp) accelerates classical Mean Shift (\#34 @comaniciu2002mean) from $O(N^2 T)$ to $O(N d 3^d)$ via hash-grid binning ($1000 times$ speedup). When global contour closure is paramount, *Normalized Cuts (Ncut)* (Alg. 104 @shi2000ncut) computes balanced spectral partitions on superpixel affinity graphs.
- *Medical Volumetric & 3D LiDAR Segmentation:* *Spatial Fuzzy C-Means (sFCM)* (Alg. 156 @chuang2006sfcm) weights fuzzy memberships $u_(i j)$ by local neighborhood aggregates $h_(i j) = sum_(r in cal(N)(x_i)) u_(r j)$, suppressing MRI Rician noise while retaining soft partial-volume tissue boundaries. For real-time 3D LiDAR after RANSAC ground removal, *DBSCAN-3D / Euclidean Cluster Extraction* (Alg. 168 @rusu2011pcl) uses voxel/k-d indexing in $bb(R)^3$ with range-adaptive radius $epsilon(r) = alpha r + epsilon_0$ (a LiDAR specialization of VDBSCAN \#50 @liu2007vdbscan) to isolate vehicles and pedestrians in $< 10$ ms.
- *Self-Supervised Semantic Image Clustering:* Deep augmentation-invariant pipelines replace pixel metrics with CNN/ViT embeddings: *SwAV* (Alg. 190 @caron2020swav) enforces cross-view prototype consistency via Sinkhorn–Knopp optimal transport; *Contrastive Clustering (CC)* (Alg. 192 @li2021contrastive) couples instance- and cluster-level InfoNCE losses @oord2018infonce; and *SCAN* (Alg. 189 @vangansbeke2020scan) decouples representation pretraining from entropy-regularized $k$-NN semantic alignment (outperforming end-to-end DeepCluster \#191 @caron2018deepcluster and IIC \#188 @ji2019iic).

#proposition(title: "SLIC O(N) Window Complexity and SwAV Anti-Collapse Entropy Guarantee", origin: "Original Proof")[
  1. *(SLIC $K$-Independence):* With superpixel stride $S = sqrt(N / K)$ and $2 S times 2 S$ local search, SLIC (Alg. 166 @achanta2012slic) inspects $O(1)$ centers per pixel, achieving exact $O(N T)$ complexity independent of $K$ (vs. $O(N K T)$ for $k$-Means \#1).
  2. *(SwAV Anti-Collapse Bound):* In SwAV (Alg. 190 @caron2020swav), solving the entropic optimal transport step over batch logits $C^top Z in bb(R)^(K times B)$,
     $ Q^* = arg max_(Q in cal(U)(B, K)) "Tr"(Q^top C^top Z) + epsilon H(Q), $
     on $cal(U)(B, K) := {Q >= 0 : Q bold(1)_B = (1 / K) bold(1)_K, Q^top bold(1)_K = (1 / B) bold(1)_B}$, forces batch-marginal prototype entropy $H(Q^* bold(1)_B) = log K$, strictly preventing trivial representation collapse.
]

#proof[
  *(1)* Lattice centers spaced by $S = sqrt(N / K)$ have spatial density $1 / S^2$; any pixel lies inside the $2 S times 2 S$ window of at most $ceil(4 S / S)^2 = 16$ centers, yielding $O(N)$ distance evaluations per iteration independent of $K$.
  *(2)* By Cuturi's Sinkhorn–Knopp theorem @cuturi2013sinkhorn, the unique maximizer $Q^* = "diag"(u) exp(C^top Z / epsilon) "diag"(v)$ satisfies row-sum feasibility $bar(q) := Q^* bold(1)_B = (1 / K) bold(1)_K$, giving exact maximum entropy
  $ H(bar(q)) = - sum_(j=1)^K (1 / K) log(1 / K) = log K, $
  ruling out degenerate single-prototype fixed points ($H(bar(q)) = 0$) at every mini-batch update.
]

#subsec("C", "Vertical 3: Natural Language Processing, Retrieval, and Large-Scale Embeddings")

*Data Invariants and Failure Modes:* Language and retrieval systems inhabit two high-dimensional regimes: (i) *sparse TF-IDF / count corpora* $X in bb(R)_(>= 0)^(n times V)$ ($V in [10^4, 10^6]$, $rho_"sp" > 0.99$) with multi-decade document-length variance; and (ii) *dense transformer embeddings* $E in bb(R)^(n times d)$ ($d in [384, 4096]$, $n in [10^6, 10^9]$) residing on $bb(S)^(d-1)$. Ambient Euclidean $k$-Means (\#1) confounds length with topic and becomes compute-bound at billion scale, while document-level *LDA* (\#83 @blei2003lda) collapses on micro-texts ($L < 15$ tokens).

*Canonical Algorithm Architecture:*
- *Directional & Hyperspherical Clustering:* *Spherical K-Means* (Alg. 12 @dhillon2001concept) maximizes sum-of-cosines on $bb(S)^(d-1)$ in $O("nnz"(X) k T)$ sparse time, while *von Mises–Fisher Mixtures (vMF)* (Alg. 89 @banerjee2005vmf) capture heterogeneous topic tightness via cluster-specific concentrations $kappa_j$.
- *Topic Modeling Spectrum:* For long multi-topic documents, *LDA* (Alg. 83 @blei2003lda) models Dirichlet simplex mixtures; *HDP* (Alg. 91 @teh2006hdp) infers $K$ nonparametrically; and *CTM* (Alg. 100 @blei2007ctm) replaces Dirichlet independence with logistic-normal topic covariance. For short texts (queries, tweets, tickets), single-topic *DMM* (Alg. 97 @yin2014dmm) and corpus-wide biterm *BTM* (Alg. 101 @yan2013btm) defeat document-level co-occurrence sparsity; *NMF* (\#92 @lee1999nmf) offers rapid additive part-based topics.
- *Billion-Scale Vector Retrieval (IVF-PQ) & Lexical Graphs:* Production ANN indexes (FAISS/ScaNN) couple coarse Voronoi partitioning via *K-Means++* (Alg. 2 @arthur2007kmeanspp), *Mini-Batch K-Means* (\#3 @sculley2010webscale), or distributed *K-Means||* (\#23 @bahmani2012scalable) ($K_"ivf" approx sqrt(n)$) with residual *Product Quantization (PQ)* running *LBG quantization* (Alg. 30 @linde1980algorithm) across $M$ orthogonal sub-spaces (compressing 1536-dim `f32` vectors $64 times$ into 96 bytes). On massive lexical/entity-resolution $k$-NN graphs without fixed $k$, *Chinese Whispers* (Alg. 111 @biemann2006chinese) executes linear $O(m T)$ weighted neighbor voting.

#proposition(title: "Spherical–vMF Equivalence and Short-Text Sparsity Collapse", origin: "Original Proof")[
  1. *(Spherical Limit):* On $bb(S)^(d-1)$, Euclidean SSE minimization (Alg. 1) is equivalent to Spherical K-Means (Alg. 12) and equals the $kappa -> infinity$ limit of equal-concentration vMF EM (Alg. 89 @banerjee2005vmf).
  2. *(Short-Text Posterior Collapse):* In symmetric $K$-topic LDA (Alg. 83) with prior $alpha bold(1)_K$ on a length-$L$ document, $|bb(E)[theta_k | z] - 1 / K| <= L / (K alpha)$, causing topic identifiability error $Omega(1 / sqrt(L))$ as $L / (K alpha) -> 0$. Pooling all $N_"bi" = sum_d binom(L_d, 2)$ corpus biterms in BTM (Alg. 101 @yan2013btm) restores $O(1 / (n L^2))$ concentration even at $L = 5$.
]

#proof[
  *(1)* On $bb(S)^(d-1)$, $norm(x_i - mu_j)_2^2 = 2 - 2 chevron.l x_i, mu_j chevron.r$, and maximizing $sum_(x in C_j) chevron.l x, mu_j chevron.r$ over $norm(mu_j)_2 = 1$ gives $mu_j^* = (sum_(x in C_j) x) / norm(sum_(x in C_j) x)_2$. In isotropic vMF EM,
  $ gamma_(i j) = (pi_j exp(kappa chevron.l x_i, mu_j chevron.r)) \
    / (sum_(r=1)^K pi_r exp(kappa chevron.l x_i, mu_r chevron.r)), $
  which converges to hard cosine assignment as $kappa -> infinity$, matching Spherical K-Means identically.
  *(2)* Conditional on topic counts summing to $L$, $theta_d | z_d tilde.op "Dir"(alpha + n_(d, 1), dots, alpha + n_(d, K))$ has mean $(alpha + n_(d, k)) / (K alpha + L)$, deviating from $1 / K$ by at most $L / (K alpha)$. Under BTM, $N_"bi" = n L (L - 1) / 2$ global pairs share one corpus proportion $theta_"corp"$, driving estimation variance down at rate $O(1 / (n L^2))$ independent of single-document brevity.
]

#subsec("D", "Vertical 4: Cybersecurity, Fraud, and Anomaly Detection")

*Data Invariants and Failure Modes:* Security telemetry (NetFlow, card transactions, endpoint logs) combines extreme class imbalance ($eta_"out" in [10^(-4), 0.05]$), massive adversarial coordinate spikes (volumetric exfiltration/DDoS), multi-density normal baselines (cron jobs vs. interactive logins), and real-time concept drift. Unprotected means and Gaussian MLEs (\#1, \#81) fail at breakdown $epsilon^* = 1 / n$ (Rule 3), creating the *masking effect* where extreme attacks inflate cluster variance and hide subsequent intrusions.

*Canonical Algorithm Architecture:*
- *Density Excision & Local Outlier Scoring:* *DBSCAN* (Alg. 31 @ester1996dbscan) and *HDBSCAN* (Alg. 32 @campello2013hdbscan) isolate anomalies via explicit noise partition $C_(-1)$ and GLOSH scores. Under multi-density network regimes, *OPTICS-OF* (Alg. 53 @breunig1999opticsof) computes local reachability outlier factors directly from the OPTICS ordering (\#33 @ankerst1999optics), flagging local density drops without global bandwidth distortion.
- *Robust Parametric Profiles:* Where regulated fraud controls require parametric behavioral baselines, *Trimmed K-Means* (Alg. 26 @cuesta1997trimmed) excises the top $alpha$-fraction of residual vectors per step, while *Student's $t$-Mixture Models (TMM)* (Alg. 88 @peel2000robust) continuously downweight extreme spikes via latent scale weights $u_(i j)$ (complemented by *PCM* \#152 @krishnapuram1993pcm and *Skew-Normal Mixtures* \#98 @lin2009skewnormal).
- *Streaming Intrusion Detection:* At $10^5$--$10^6$ events/sec, *DenStream* (Alg. 194 @cao2006denstream) separates exponentially damped core ($p$) and outlier ($o$) micro-clusters—detecting emerging botnet waves as $o$-clusters promote across threshold $beta mu$; *D-Stream* (Alg. 197 @chen2007dstream) isolates sporadic grid cells in $O(1)$ online time; and *CluStream* (Alg. 193 @aggarwal2003clustream) enables pyramidal forensic time-horizon comparisons.

#theorem(title: "Redescending Influence Bound of Student-t Mixtures and Trimmed Masking Immunity", origin: "Original Proof")[
  Let $cal(X) = cal(X)_"clean" union cal(X)_"adv"$ contain $(1 - eta) n$ clean points and $eta n$ arbitrary adversarial outliers ($eta < 1 / 2$).
  1. *(TMM Redescending Leverage):* In TMM (Alg. 88 @peel2000robust) with degrees of freedom $nu_j > 0$ and $Delta_(i j)^2 := (x_i - mu_j)^top Sigma_j^(-1) (x_i - mu_j)$, the centroid influence $psi_j (x_i) := gamma_(i j) u_(i j) (x_i - mu_j)$ ($u_(i j) = (nu_j + d) / (nu_j + Delta_(i j)^2)$) is uniformly bounded and strictly redescending:
     $ norm(Sigma_j^(-1 / 2) psi_j (x_i))_2 <= (nu_j + d) / (2 sqrt(nu_j)), \
       lim_(Delta_(i j) -> infinity) norm(Sigma_j^(-1 / 2) psi_j (x_i))_2 = 0, $
     whereas GMM (\#81) suffers unbounded influence $gamma_(i j) Delta_(i j) -> infinity$.
  2. *(Trimmed Zero-Masking Guarantee):* With $alpha >= eta$ initialized inside the clean basin of radius $R_"clean"$, every outlier with $min_j norm(x_"adv" - mu_j)_2 > R_"clean"$ exerts identically zero gradient on all $k$ centroids of Trimmed K-Means (Alg. 26 @cuesta1997trimmed).
]

#proof[
  *(1)* With $gamma_(i j) in [0, 1]$ and whitened norm $norm(Sigma_j^(-1 / 2) (x_i - mu_j))_2 = Delta_(i j)$,
  $ norm(Sigma_j^(-1 / 2) psi_j (x_i))_2 &= gamma_(i j) ((nu_j + d) Delta_(i j)) / (nu_j + Delta_(i j)^2) \
    &<= ((nu_j + d) Delta_(i j)) / (nu_j + Delta_(i j)^2). $
  AM-GM gives $nu_j + Delta_(i j)^2 >= 2 sqrt(nu_j) Delta_(i j)$ (tight at $Delta_(i j) = sqrt(nu_j)$), yielding ceiling $(nu_j + d) / (2 sqrt(nu_j))$, while dividing by $Delta_(i j)^2$ drives the limit to $0$ as $Delta_(i j) -> infinity$. Covariance leverage $u_(i j) Delta_(i j)^2 <= nu_j + d$ also remains uniformly bounded.
  *(2)* Since $|cal(X)_"clean"| = (1 - eta) n >= (1 - alpha) n$ and every clean squared residual is $<= R_"clean"^2 < d_"adv"$, the C-step order statistic $d_((floor((1 - alpha) n))) <= R_"clean"^2$ strictly excludes $cal(X)_"adv"$ from the active set $A_alpha$, yielding $partial hat(mu)_j / partial x_"adv" = 0$.
]

#subsec("E", "Vertical 5: Quantitative Finance, Portfolio Risk, and Customer Analytics")

*Data Invariants and Failure Modes:* Quantitative finance and CRM span three modalities: (i) *mixed numeric–categorical tables* pairing heavy-tailed financial ratios with nominal industry/product codes; (ii) *cross-sectional asset correlation matrices* ($n in [50, 5000]$, $d in [10, 250]$) with elliptical factor collinearity and fat tails; and (iii) *price, flow, and spend time series* exhibiting lead–lag phase shifts and volatility warping. Raw $L_2$ on one-hot CRM codes lets high-cardinality nominals swamp continuous risk variables, while pointwise $L_2$ misclassifies identical sector peers shifted by a 2-day earnings lag.

*Canonical Algorithm Architecture:*
- *Mixed CRM & Credit Segmentation:* *K-Prototypes* (Alg. 8 @huang1998extensions) combines numeric $L_2$ and categorical Hamming cost ($gamma$) in $O(n k d T)$ time. When arbitrary business dissimilarities (Gower, Jaccard) and *real exemplar medoids* are required, *PAM* (Alg. 4 @kaufman1990finding) and its $O(k)$-accelerated exact successor *FasterPAM* (Alg. 24 @schubert2019fasterpam) deliver robust, auditable prototypes. Under noisy marketing attributes, *W-K-Means* (Alg. 28 @huang2005automated) and *Minkowski Weighted K-Means* (\#29 @deamorim2012minkowski) learn convex feature weights automatically; probabilistic alternatives include *LCA* (\#90 @lazarsfeld1968latent) and *Mixmod* (\#94 @biernacki2006mixmod).
- *Portfolio Risk & Regime Discovery:* *Ward's Method* (Alg. 59 @ward1963hierarchical) on correlation distances $d_(i j) = sqrt((1 - rho_(i j)) / 2)$ drives *Hierarchical Risk Parity (HRP)* allocation, avoiding unstable Markowitz $Sigma^(-1)$ inversion. For market regimes, *MCLUST* (Alg. 93 @fraley2002mclust) selects among 14 eigen-decomposed covariance geometries via BIC @schwarz1978bic, balancing correlation flexibility against finite-sample overfitting.
- *Financial Time-Series Regimes:* *k-Shape* (Alg. 162 @paparrizos2015kshape) computes shift- and scale-invariant *Shape-Based Distance (SBD)* via $O(m log m)$ FFT cross-correlation. Under nonlinear time warping, *Soft-DTW Clustering* (Alg. 164 @cuturi2017softdtw) replaces non-smooth *DBA* (\#163 @petitjean2011dba) with differentiable $min_gamma$ barycenters; *ROCKET-Clustering* (\#165 @dempster2020rocket) provides a linear-time random-kernel baseline for long tick streams.

#proposition(title: "MCLUST Eigen-Covariance Parsimony and k-Shape Shift-Scale Invariance", origin: "Original Proof")[
  1. *(MCLUST Parameter Economy):* Eigen-decomposing $Sigma_j = lambda_j D_j A_j D_j^top$ in MCLUST (Alg. 93 @fraley2002mclust) compresses covariance parameters from $p_"VVV" = k d (d + 1) / 2$ to $p_"EEE" = d (d + 1) / 2$, $p_"VVI" = k d$, or $p_"EII" = 1$. Whenever $n < 5 k d^2 / 2$, BIC rejects unconstrained VVV unless genuine cross-regime eigen-rotations supply $Omega(d log n / n)$ average log-likelihood gain per sample.
  2. *(Shift-Affine Invariance of SBD & Smoothness of Soft-DTW):* For non-constant $x, y in bb(R)^m$, $a_1, a_2 > 0$, $b_1, b_2 in bb(R)$, and circular shift $S_s$, k-Shape (Alg. 162 @paparrizos2015kshape) satisfies
     $ "SBD"(a_1 S_s (x) + b_1 bold(1), a_2 y + b_2 bold(1)) = "SBD"(x, y) in [0, 2], $
     while Soft-DTW (Alg. 164 @cuturi2017softdtw, $gamma > 0$) is $C^infinity$-smooth with gradient $bb(E)_(P_gamma) [A] (x - y)$.
]

#proof[
  *(1)* Orthogonal $D_j in O(d)$ holds $d (d - 1) / 2$ parameters, unit-determinant diagonal $A_j$ holds $d - 1$, and $lambda_j$ holds $1$, giving $d (d + 1) / 2$ per component under VVV vs. $1$ shared block under EEE—saving $(k - 1) d (d + 1) / 2$ parameters. Under a shared-covariance null, Wilks' likelihood ratio scales as $O_p (k d^2)$ while BIC penalizes VVV by $((k - 1) d (d + 1) / 4) log n$, preventing covariance overfitting unless per-sample rotation evidence exceeds $Omega(d log n / n)$.
  *(2)* Z-normalization $tilde(x) = (x - bar(x) bold(1)) / norm(x - bar(x) bold(1))_2$ eliminates $a_1 > 0$ and $b_1$ identically; maximizing normalized cross-correlation $"CC"_w$ over all lags $w$ absorbs $S_s$ by re-indexing $w |-> w - s$, with $|"CC"_w| <= 1$ by Cauchy–Schwarz. Differentiating $"dtw"_gamma = -gamma log sum_A exp(-chevron.l A, Delta chevron.r / gamma)$ yields $nabla_Delta "dtw"_gamma = bb(E)_(P_gamma) [A]$ smoothly for all $gamma > 0$.
]

#subsec("F", "Vertical 6: High-Velocity IoT, Network Telemetry, and Edge Streaming")

*Data Invariants and Failure Modes:* Industrial IoT meshes, 5G packet brokers, and Kafka/Flink edge buses emit unbounded streams ${x_t}_(t=1)^infinity subset.eq bb(R)^d$ ($10^4$--$10^7$ events/sec) under strict single-pass arrival, bounded gateway RAM ($B_"mem" = o(n d)$), erratic inter-arrival bursts ("anytime" interruption), and continuous concept drift. Every multi-pass batch method violates Rule 1 of Theorem 1.

*Canonical Algorithm Architecture:*
- *Additive CF-Trees & Coreset Summarizers:* *BIRCH* (Alg. 64 @zhang1996birch) compresses massive Euclidean data in one $O(n d)$ pass via a height-balanced CF-tree storing exact additive triplets $"CF" = (N, bold("LS"), "SS")$ within fixed RAM. *BICO* (Alg. 78 @fichtenberger2013bico) replaces BIRCH's heuristic radius with *provable $(1 + epsilon)$-coreset guarantees* via nearest-neighbor projection, while *StreamKM++* (Alg. 195 @ackermann2012streamkmpp) applies Bentley--Saxe *merge-and-reduce* over K-Means++ coreset trees (\#2) in polylogarithmic space.
- *Online–Offline & Anytime Evolving Streams:* *CluStream* (Alg. 193 @aggarwal2003clustream) maintains online timestamped micro-clusters and pyramidal snapshots enabling offline queries over arbitrary historical windows $[t_1, t_2]$. *ClusTree* (Alg. 196 @kranen2011clustree) supports *anytime* insertion via hitchhiking R-tree descent with embedded exponential decay, dropping zero packets during burst floods.
- *Non-Convex Density/Grid Streams & Mini-Batch SGD:* *DenStream* (Alg. 194 @cao2006denstream) maintains damped core and outlier micro-clusters for arbitrary-shape drift, while *D-Stream* (Alg. 197 @chen2007dstream) updates decaying density grids in $O(1)$ online time ($d <= 12$). For spherical centroids at maximum throughput, *Mini-Batch K-Means* (Alg. 3 @sculley2010webscale) executes $O(b k d)$ stochastic gradient steps @bottou1995convergence.

#theorem(title: "CF Additivity, Fading Memory Ceiling, and Merge-and-Reduce Coreset Bound", origin: "Original Unified Proof")[
  Let ${x_t}_(t=1)^T subset.eq bb(R)^d$ arrive sequentially with exponential decay $2^(-lambda Delta t)$ ($lambda > 0$, unit arrival rate).
  1. *(CF Sufficiency — BIRCH \#64 / CluStream \#193 / ClusTree \#196):* For disjoint sets $A, B$, $"CF"(A union B) = "CF"(A) + "CF"(B)$ with $"CF"(A) := (N_A, bold("LS")_A, "SS"_A)$, yielding exact radius and SSE in $O(d)$ time without raw points:
     $ R(A)^2 = "SS"_A / N_A - norm(bold("LS")_A / N_A)_2^2, \
       "SSE"(A) = "SS"_A - norm(bold("LS")_A)_2^2 / N_A. $
  2. *(Constant-RAM Ceiling — DenStream \#194 / D-Stream \#197):* Under exponential damping, total system weight is bounded uniformly by $W_"max" = 1 / (1 - 2^(-lambda))$, capping active core micro-clusters (weight $>= beta mu$) across an infinite stream at @cao2006denstream, @chen2007dstream
     $ M_"max" <= floor( 1 / (beta mu (1 - 2^(-lambda))) ) = O(1). $
  3. *(StreamKM++ \#195 / BICO \#78 Bound):* Using $L = ceil(log_2 (T / m))$ merge-and-reduce buffers of size $m = O(k d epsilon^(-2) log(k / delta) log^2 T)$ plus K-Means++ (\#2) guarantees expected cost $<= (1 + epsilon) dot 8 (ln k + 2) "OPT"_k (cal(X)_T)$ @ackermann2012streamkmpp, @arthur2007kmeanspp.
]

#proof[
  *(1)* Expanding $"SSE"(A) = sum_(x in A) norm(x - mu_A)_2^2$ with $mu_A = bold("LS")_A / N_A$ gives
  $ "SSE"(A) &= "SS"_A - 2 chevron.l bold("LS")_A, bold("LS")_A / N_A chevron.r \
    &quad + N_A norm(bold("LS")_A / N_A)_2^2 = "SS"_A - norm(bold("LS")_A)_2^2 / N_A, $
  and linearity of summation yields exact additivity across disjoint sets.
  *(2)* At time $t$, unit arrivals accumulate total damped weight
  $ W(t) &= sum_(s=0)^(t-1) 2^(-lambda s) \
    &= (1 - 2^(-lambda t)) / (1 - 2^(-lambda)) < 1 / (1 - 2^(-lambda)) = W_"max". $
  Since disjoint core clusters each hold mass $>= beta mu$, $M(t) dot beta mu <= W(t) < W_"max"$, bounding $M(t)$ uniformly for all $t$.
  *(3)* By coreset composability @harpeled2004coresets, setting per-level error $epsilon' = epsilon / (2 L)$ across depth $L$ incurs total factor $(1 + epsilon / (2 L))^L <= exp(epsilon / 2) <= 1 + epsilon$; composing with Arthur--Vassilvitskii's $8 (ln k + 2)$ K-Means++ bound @arthur2007kmeanspp completes the proof.
]

#subsec("G", "Vertical 7: Geospatial Analytics, Urban Mobility, and Spatiotemporal Tracking")

*Data Invariants and Failure Modes:* GIS grids, connected-vehicle GPS pings, and maritime AIS tracks combine (i) *heterogeneous product metrics* coupling geodesic coordinates $(phi, lambda)$, timestamp $t$, and continuous attributes (temperature, speed); (ii) *orders-of-magnitude urban–rural density gradients*; and (iii) *variable-length polylines* sharing common corridor sub-segments. Naive $L_2$ concatenation of $(x, y, t)$ fails unit-scale invariance, while fixed-$epsilon$ DBSCAN (\#31) either merges entire metropolitan areas or erases rural corridors.

*Canonical Algorithm Architecture:*
- *Spatiotemporal & Multi-Density Point Extraction:* *ST-DBSCAN* (Alg. 40 @birant2007stdbscan) decouples spatial radius $epsilon_1$, temporal tolerance $epsilon_2$, and cluster-mean drift guard $Delta epsilon$ to prevent seasonal chaining. *OPTICS* (Alg. 33 @ankerst1999optics) exposes nested urban-to-street reachability hierarchies, while *VDBSCAN* (Alg. 50 @liu2007vdbscan) extracts multi-strata $k$-dist knees $epsilon_1 < dots < epsilon_m$ sequentially. *GDBSCAN* (\#39 @sander1998gdbscan) handles polygon/road-network predicates, and *GriDBSCAN* (\#48 @mahran2008gridbscan) / *MR-DBSCAN* (\#49 @he2011mrdbscan) scale across distributed Spark clusters.
- *Trajectory Sub-Segment & Multiresolution GIS Engines:* *TRACLUS* (Alg. 169 @lee2007traclus) partitions GPS tracks at MDL turning points and density-clusters directed line segments under perpendicular, parallel, and angular distances ($d_perp + d_parallel + d_theta$). *STING* (Alg. 133 @wang1997sting) builds statistical quadtrees answering spatial queries in $O(K_"cells")$ time independent of $n$; *WaveCluster* (Alg. 135 @sheikholeslami1998wavecluster) applies multi-resolution discrete wavelets to spatial grids in $O(n)$ time, simultaneously denoising GPS jitter; and *CLARANS* (Alg. 6 @ng2002clarans) scales spatial $k$-medoid hub siting via randomized swap-graph search.

#proposition(title: "Product-Metric Decoupling in ST-DBSCAN and MDL Optimality of TRACLUS", origin: "Original Proof")[
  1. *(Unit-Scale Invariance of ST-DBSCAN \#40):* On $p = (s_p, t_p, a_p)$, concatenated $L_2$ balls fail under independent physical unit changes $(s, t, a) |-> (c_s s, c_t t, c_a a)$ ($c_s != c_t$), whereas ST-DBSCAN's product predicate @birant2007stdbscan
     $ N_(epsilon_1, epsilon_2) (p) := { & q : d_"geo" (s_p, s_q) <= epsilon_1 \
       & and |t_p - t_q| <= epsilon_2 } $
     with running mean guard $|a_q - bar(a)_C| <= Delta epsilon$ is strictly equivariant under diagonal unit scaling and blocks monotonic attribute chaining across length $L > Delta epsilon / Delta_"step"$.
  2. *(TRACLUS MDL Criterion \#169):* Merging polyline $p_1, dots, p_m$ into macro-segment $L_(1 m)$ compresses description length iff $L(H) + L(D | H) < L_"nopar"$ @lee2007traclus, balancing segment compactness against perpendicular and angular turning fidelity.
]

#proof[
  *(1)* Rescaling meters to kilometers ($c_s = 10^(-3), c_t = 1$) collapses spatial contribution by $10^6$ in concatenated $L_2$, reducing clustering to pure 1D time slicing. In ST-DBSCAN, each factor inequality is degree-1 homogeneous in its own physical dimension, ensuring exact equivariance under $(epsilon_1, epsilon_2, Delta epsilon) |-> (c_s epsilon_1, c_t epsilon_2, c_a Delta epsilon)$; meanwhile $|a_q - bar(a)_C| <= Delta epsilon$ halts linear drift $a_r = r delta$ once deviation from the cluster mean exceeds $Delta epsilon$.
  *(2)* Raw prefix coding requires $L_"nopar" = sum_(j=1)^(m-1) log_2 norm(p_(j+1) - p_j)_2$ bits, whereas macro-segment coding costs $L(H) = log_2 norm(p_m - p_1)_2$ plus residual deviation bits $L(D | H) = log_2 (sum d_perp) + log_2 (sum d_theta)$. Collinear paths satisfy $L(D | H) = 0$ and merge greedily, while sharp turns inflate $L(D | H)$ until $L(H) + L(D | H) > L_"nopar"$, forcing an exact characteristic split @lee2007traclus.
]

#subsec("H", "Vertical 8: Social, Web, Communication, and Biological Network Science")

*Data Invariants and Failure Modes:* Social, web, telecom, and protein-interaction graphs $G = (V, E, W)$ ($|V| in [10^3, 10^9]$) exhibit sparse power-law degree distributions ($P(k) prop k^(-gamma)$), pervasive overlapping affiliations, directed flow asymmetry, and occasional bipartite/core-periphery disassortative blocks. Modularity maximization suffers the *resolution limit* on large graphs, while the standard *SBM* (\#130) fatally clusters vertices by degree rather than community on heavy-tailed hub networks.

*Canonical Algorithm Architecture:*
- *Large-Scale Disjoint & Flow Communities:* *Leiden* (Alg. 108 @traag2019leiden) supersedes *Louvain* (Alg. 107 @blondel2008louvain) in near-linear $O(m)$ time with connected-community guarantees and Constant Potts Model (CPM) resolution control. For directed diffusion networks (citations, clickstreams), *Infomap* (Alg. 113 @rosvall2008infomap) minimizes random-walk Map Equation coding length. On biological PPI/homology graphs, *MCL* (Alg. 109 @vandongen2000mcl) alternates stochastic expansion and inflation to recover tight protein complexes.
- *Generative Blockmodels, Overlap & Balanced Cuts:* *DC-SBM* (Alg. 131 @karrer2011dcsbm) repairs standard *SBM* (\#130 @holland1983sbm) on scale-free graphs and captures disassortative/core-periphery roles. For overlapping communities, *BigCLAM* (Alg. 125 @yang2013bigclam) scales non-negative affiliation factorization in $O(m)$ per sweep (capturing dense overlaps), *OSLOM* (Alg. 124 @lancichinetti2011oslom) optimizes configuration-null $p$-values while flagging homeless noise nodes, and *CPM* (Alg. 127 @palla2005cpm) rolls adjacent $k$-cliques on biological networks (alongside *SLPA* \#126 @xie2011slpa and *DEMON* \#128 @coscia2012demon). For strictly balanced distributed database sharding and VLSI meshes, *METIS* (Alg. 121 @karypis1998metis) executes $O(m)$ multilevel coarsening with Kernighan–Lin refinement (\#120).

#theorem(title: "Modularity Resolution Limit, DC-SBM Degree Separation, and Kesten–Stigum Threshold", origin: "Classical — Cited")[
  Let $G = (V, E)$ be an unweighted graph with $m = |E|$ edges.
  1. *(Resolution Limit @fortunato2007resolution):* Under Newman–Girvan modularity (Louvain \#107, modularity-Leiden \#108), merging two distinct $l$-edge cliques joined by one bridge strictly increases $Q$ whenever $l < sqrt(m / 2)$.
  2. *(DC-SBM Degree Decoupling @karrer2011dcsbm):* Standard SBM (\#130) groups hubs apart from leaves whenever degree heterogeneity dominates block contrast; DC-SBM (Alg. 131) yields profile log-likelihood
     $ cal(L)_"DC-SBM" (g) = sum_(r, s) m_(r s) log( m_(r s) / (kappa_r kappa_s) ), $
     with block volume $kappa_r = sum_(i in C_r) k_i$, decoupling vertex degree completely from community affinity.
  3. *(Kesten–Stigum Detectability @abbe2018sbmreview):* In a symmetric 2-block sparse SBM with edge probabilities $a / n$ (intra) and $b / n$ (inter), weak recovery above chance is possible iff
     $ (a - b)^2 > 2 (a + b). $
]

#proof[
  *(1)* Merging two $l$-edge cliques (each of total degree $2 l + 1$) joined by one edge changes modularity by
  $ Delta Q &= (2 l + 1) / m - ((4 l + 2) / (2 m))^2 \
    &quad - 2 ( l / m - ((2 l + 1) / (2 m))^2 ) \
    &= 1 / m - ((2 l + 1)^2) / (2 m^2), $
  which is strictly positive iff $(2 l + 1)^2 < 2 m$, i.e., $l < (sqrt(2 m) - 1) / 2 approx sqrt(m / 2)$ @fortunato2007resolution.
  *(2)* Maximizing the Poisson DC-SBM log-likelihood
  $ cal(L)(theta, omega, g) &= sum_(i, j) A_(i j) log(theta_i theta_j omega_(g_i g_j)) \
    &quad - sum_(i, j) theta_i theta_j omega_(g_i g_j) $
  subject to $sum_(i in C_r) theta_i = 1$ yields exact MLEs $hat(theta)_i = k_i / kappa_(g_i)$ and $hat(omega)_(r s) = m_(r s)$; substituting back cancels degree factors into constants, leaving $sum_(r, s) m_(r s) log(m_(r s) / (kappa_r kappa_s))$ @karrer2011dcsbm. Normalizing by volume $kappa_r kappa_s$ removes hub–leaf likelihood bias.
  *(3)* The non-backtracking operator has informative eigenvalue $lambda_2 = (a - b) / 2$ against bulk spectral radius $sqrt((a + b) / 2)$; escaping the semicircle requires $(a - b) / 2 > sqrt((a + b) / 2) <=> (a - b)^2 > 2 (a + b)$ @abbe2018sbmreview.
]

#subsec("I", "Vertical 9: High-Dimensional Industrial Sensors and Subspace Discovery")

*Data Invariants and Failure Modes:* Semiconductor wafer-fab telemetry ($d in [200, 5000]$ sensors), aerospace turbines, and hyperspectral arrays exhibit *cluster-specific low-dimensional subspaces*: each fault mode lives inside $d_j in [3, 15]$ relevant or coupled sensors while all remaining coordinates emit process noise. Global PCA mixes distinct fault subspaces and buries low-variance fault axes beneath background drift, while ambient $L_2$ fails by Rule 2 of Theorem 1.

*Canonical Algorithm Architecture:*
- *Axis-Parallel & Statistical Subspace Discovery:* For strict physical-sensor interpretability ($iota_"int" = "strict"$), *CLIQUE* (Alg. 134 @agrawal1998clique) and *SUBCLU* (Alg. 137 @kailing2004subclu) exhaustively enumerate dense axis-aligned subspaces via Apriori pruning (practical for $d_"max" <= 12$, Chapter XII Theorem 1). To eliminate exponential lattice explosion and redundant sub-projections, *P3C* (Alg. 148 @moise2006p3c) merges statistically significant 1D Poisson intervals into cluster cores, while *STATPC* (Alg. 149 @moise2008statpc) extracts non-redundant statistically significant clusters against a uniform null (alongside *FIRES* \#146 @kriegel2005fires and *INSCY* \#147 @assent2008inscy).
- *Top-Down Projected & Soft Attribute Weighting:* *PROCLUS* (Alg. 136 @aggarwal1999proclus) assigns every observation to one disjoint medoid-centered axis subspace of target dimension $ell$ in $O(n k d)$ time; *COSA* (Alg. 145 @friedman2004cosa) learns smooth exponential-entropy coordinate weights per cluster without hard attribute dropout.
- *Oblique Linear Subspaces:* When physical conservation laws couple sensors along rotated manifolds, *ORCLUS* (Alg. 138 @aggarwal2000orclus) discovers arbitrarily oriented subspaces via local cluster PCA, while *Sparse Subspace Clustering (SSC)* (Alg. 150 @elhamifar2013ssc) solves the convex $ell_1$ self-expression program
  $ min_(c_i) norm(c_i)_1 quad "s.t." quad x_i = X c_i, quad c_(i i) = 0, $
  and partitions $W = |C| + |C|^top$ spectrally—provably separating intersecting, rotated subspaces of unknown dimensions.

#proposition(title: "SSC Subspace-Detection Exactness and Axis-Rotation Blindness", origin: "Classical — Cited")[
  1. *(Axis-Parallel Obstruction):* Two orthogonal 1D clusters along $u_1 = d^(-1 / 2) bold(1)$ and $u_2 = d^(-1 / 2) (1, -1, 1, -1, dots)^top$ with equal line variance $sigma^2$ share *identical* univariate marginals $cal(N)(0, sigma^2 / d)$ on every coordinate $r = 1, dots, d$, blinding every 1D axis-parallel method (CLIQUE \#134, PROCLUS \#136, P3C \#148).
  2. *(SSC Geometric Exactness @elhamifar2013ssc, @soltanolkotabi2012geometric):* For unit-norm samples from subspaces $cal(U)_1, dots, cal(U)_k$, if every point's in-subspace convex-hull inradius $r(cal(P)_(-i)^((l)))$ exceeds maximum dual coherence to external subspaces,
     $ max_(j in.not cal(X)_l) |chevron.l v_i^((l)), x_j chevron.r| < r(cal(P)_(-i)^((l))), $
     then every SSC optimum $c_i^*$ is strictly subspace-preserving ($c_(i j)^* = 0$ for $x_j in.not cal(U)_l$).
]

#proof[
  *(1)* Both $z u_1$ and $w u_2$ ($z, w tilde.op cal(N)(0, sigma^2)$) project onto coordinate $r$ as $plus.minus z / sqrt(d) tilde.op cal(N)(0, sigma^2 / d)$, erasing all 1D marginal contrast despite $90 degree$ subspace orthogonality, whereas local PCA (ORCLUS \#138) and $ell_1$ self-expression (SSC \#150) exploit cross-coordinate correlation directly.
  *(2)* Restricting SSC to $cal(U)_l$ yields primal optimum $c_i^((l))$ and dual certificate $nu_i in cal(U)_l$ with $norm(nu_i)_2 <= 1 / r(cal(P)_(-i)^((l)))$ by polarity of $cal(P)_(-i)^((l)) = "conv"(plus.minus cal(X)_(-i)^((l)))$. Setting $v_i^((l)) := nu_i / norm(nu_i)_2$, for every external $x_j in.not cal(X)_l$:
  $ |chevron.l nu_i, x_j chevron.r| &= norm(nu_i)_2 |chevron.l v_i^((l)), x_j chevron.r| \
    &<= |chevron.l v_i^((l)), x_j chevron.r| / r(cal(P)_(-i)^((l))) < 1, $
  satisfying strict Basis Pursuit KKT optimality with zero inter-subspace edges @soltanolkotabi2012geometric.
]

#subsec("J", "Vertical 10: Multi-Modal Fusion and Unsupervised Consensus Pipelines")

*Data Invariants and Failure Modes:* Multi-omics medicine (RNA + methylation + pathology), multimodal commerce (text + image + co-purchase graph), and sensor fusion observe $n$ entities across $V >= 2$ heterogeneous views ${X^((1)), dots, X^((V))}$. Raw feature concatenation fails via dimensionality dominance, single-view noise contamination, and metric incompatibility.

*Canonical Algorithm Architecture:*
- *Joint Multi-View Models:* *Multi-View K-Means* (Alg. 174 @cai2013multiview) learns adaptive view weights $alpha_v^gamma$ under robust $ell_(2, 1)$ distortion, automatically attenuating noisy views. *Co-regularized Multi-View Spectral Clustering* (Alg. 175 @kumar2011coregularized) aligns view-specific Laplacian eigenspaces via Hilbert–Schmidt regularization to fuse non-convex manifold geometries.
- *Label-Only Ensemble Consensus:* Given $R$ diverse partitions (across views or families), Strehl and Ghosh @strehl2002cluster fuse labelings without raw features via *CSPA* (Alg. 171, $O(n^2 k R)$ co-association METIS cut for $n <= 10^4$), *HGPA* (Alg. 172, $O(n k R)$ hypergraph cut), and *MCLA* (Alg. 173, $O(n k^2 R^2)$ meta-cluster Jaccard partitioning scaling cleanly to millions of points).
- *Deep Latent Fusion:* *VaDE* (Alg. 184 @jiang2017vade) embeds a Gaussian Mixture prior directly inside a VAE ELBO @kingma2014vae for calibrated multimodal uncertainty, while *IDEC* (Alg. 182 @guo2017idec) retains autoencoder reconstruction alongside DEC's KL self-training (\#181 @xie2016dec) to preserve local latent structure.

#theorem(title: "Consensus Variance Contraction and Co-Regularized Spectral Perturbation Bound", origin: "Original Unified Proof")[
  1. *(Co-Association Variance Contraction — Monti \#170 / CSPA \#171):* For $R$ base partitions with binary pair indicators $A_(i j)^((r))$, variance $v_(i j) = p_(i j)(1 - p_(i j))$, and average inter-run correlation $bar(rho)_(i j) in [0, 1)$, the consensus matrix $bar(S)_(i j) := (1 / R) sum_(r=1)^R A_(i j)^((r))$ satisfies
     $ "Var"(bar(S)_(i j)) = v_(i j) ( bar(rho)_(i j) + (1 - bar(rho)_(i j)) / R ), $
     minimized by maximizing cross-paradigm diversity (lowering $bar(rho)_(i j)$).
  2. *(Multi-View Spectral Contraction — Alg. 175 @kumar2011coregularized):* If $V$ views observe shared $k$-cluster Laplacian $L_*$ (eigengap $delta > 0$) via $hat(L)^((v)) = L_* + E^((v))$ with independent zero-mean noise $norm(E^((v)))_2 <= sigma_v$, centroid co-regularization ($lambda_v -> infinity$) achieves Davis–Kahan subspace error @davis1970rotation
     $ bb(E)[norm(sin Theta(hat(U)^*, U_*))_F] \
       &quad <= (sqrt(2 k)) / (delta V) sqrt(sum_(v=1)^V bb(E)[norm(E^((v)))_2^2]) \
       &quad = O( bar(sigma) / (delta sqrt(V)) ), $
     contracting spectral recovery error strictly by $sqrt(V)$.
]

#proof[
  *(1)* Expanding the variance of $bar(S)_(i j)$ gives
  $ "Var"(bar(S)_(i j)) &= 1 / R^2 ( sum_(r=1)^R "Var"(A_(i j)^((r))) \
    &quad + sum_(r != s) "Cov"(A_(i j)^((r)), A_(i j)^((s))) ) \
    &= v_(i j) ( bar(rho)_(i j) + (1 - bar(rho)_(i j)) / R ). $
  Multi-seed restarts of one biased algorithm exhibit $bar(rho)_(i j) approx 0.9$, hitting an irreducible variance floor, whereas orthogonal families suppress $bar(rho)_(i j)$ on spurious splits while keeping $p_(i j) approx 1$ on true clusters.
  *(2)* In centroid co-regularized spectral clustering minimizing
  $ & sum_(v=1)^V "Tr"(U^((v) top) hat(L)^((v)) U^((v))) \
    &quad - sum_(v=1)^V lambda_v "Tr"(U^((v)) U^((v) top) U^* U^(* top)), $
  sending $lambda_v -> infinity$ enforces common projection $U^* U^(* top)$, reducing the problem to eigen-decomposition of $bar(L) = L_* + bar(E)$ ($bar(E) = (1 / V) sum_v E^((v))$). By the Davis–Kahan $sin Theta$ theorem @davis1970rotation, $norm(sin Theta(hat(U)^*, U_*))_F <= (sqrt(2 k) / delta) norm(bar(E))_2$; independence and zero mean of $E^((v))$ give $bb(E)[norm(bar(E))_2^2] = O(V^(-2) sum_v bb(E)[norm(E^((v)))_2^2])$, and Jensen's inequality yields the $O(bar(sigma) / (delta sqrt(V)))$ bound.
]

// =============================================================================
// SECTION III: MASTER DIAGNOSTIC DECISION MATRICES
// =============================================================================
#sec("III", "Master Diagnostic Decision Matrices Across Verticals and Regimes")

Tables I and II distill the ten vertical engineering architectures of Section II and the six elimination rules of Theorem 1 into two orthogonal lookup references: *by application vertical* (Table I) and *by observable structural data invariant* (Table II). Notably, fewer than 35 core workhorse algorithms from $cal(A)_(200)$ cover over 90% of Tier-1 production workloads, while the remaining 165 methods serve as high-value domain specialists (e.g., NJ \#75, TRACLUS \#169, METIS \#121, Bimax \#159) or algebraic accelerators (FasterPAM \#24, Mean Shift++ \#55, Elkan/Hamerly \#13--14).

#wide_figure(caption: [*Table I:* Master Ten-Vertical Domain-to-Algorithm Engineering Matrix. Maps each vertical directly to exact catalog IDs from Parts I--VIII, contrasting Tier-1 primary candidates against scalable/specialized fallbacks, dominant naive failure modes, and domain validation protocols.])[
  #set par(justify: false, first-line-indent: 0pt)
  #set text(size: 6.7pt)
  #show math.equation: set text(size: 6.7pt)
  #table(
    columns: (1.18fr, 0.95fr, 1.15fr, 1.35fr, 1.12fr, 1.15fr),
    align: (left, left, left, left, left, left),
    inset: (x: 3.0pt, y: 2.5pt),
    fill: (_, y) => if y == 0 { rgb("#f1f5f9") } else if calc.odd(y) { rgb("#fbfdff") } else { white },
    table.header(
      [Scientific / Industrial Vertical], [Typical Data Invariants], [Why Naive Baselines Fail], [Primary Tier-1 Candidates (Exact IDs)], [Scalable / Niche Fallbacks (IDs)], [Critical Tuning & Validation Protocol],
    ),
    [*1. Genomics & Systems Biology*], [$n in [10^2, 10^7]$, $d approx 2 times 10^4$, $rho_"sp" > 0.85$, counts / trees], [Dropout noise, library-size norm, rare cell types, rate heterogeneity], [Leiden (\#108), SNN (\#36), HDBSCAN (\#32), NJ (\#75), FABIA (\#161), Monti (\#170)], [Louvain (\#107), UPGMA (\#58), Cheng–Church (\#157), Kluger (\#123), ISA (\#160)], [Log1p/HVG/PCA($30$); marker enrichment; Monti CDF stability; CPM $gamma$],
    [*2. Vision, Medical Imaging & 3D LiDAR*], [2D/3D grids ($10^5$--$10^8$), LiDAR $bb(R)^3$ ($<20$ ms), RGB manifolds], [Pixel salt-and-pepper noise, $r^(-2)$ LiDAR sparsity, semantic invariance], [SLIC (\#166), sFCM (\#156), DBSCAN-3D (\#168), Mean Shift++ (\#55), SwAV (\#190), SCAN (\#189)], [Watershed (\#167), Ncut (\#104), CC (\#192), Quickshift (\#41), VDBSCAN (\#50)], [Spatial weight $m$ (SLIC); range-scaled $epsilon(r)$ (LiDAR); Dice/IoU & Sinkhorn entropy],
    [*3. NLP, Retrieval & Embeddings*], [Sparse TF-IDF ($V <= 10^6$, $rho_"sp" > 0.99$); dense $bb(S)^(d-1)$ ($n <= 10^9$)], [Length-norm confounding, short-text $L < 15$ LDA collapse, $O(n k d)$ cost], [Spherical KM (\#12), vMF (\#89), LDA (\#83), CTM (\#100), DMM (\#97), BTM (\#101), IVF-PQ (\#2/\#30)], [HDP (\#91), NMF (\#92), Mini-Batch (\#3), K-Means|| (\#23), Chinese Whispers (\#111)], [$L_2$-normalize; recall\@10 vs. compression (IVF-PQ); NPMI coherence; perplexity],
    [*4. Cybersecurity, Fraud & Anomaly*], [Streaming $10^5$/s, $eta_"out" in [10^(-4), 0.05]$, multi-density drift], [$1 / n$ mean breakdown, outlier masking effect, adversarial drift], [HDBSCAN (\#32), OPTICS-OF (\#53), DenStream (\#194), TMM (\#88), Trimmed KM (\#26)], [DBSCAN (\#31), CluStream (\#193), D-Stream (\#197), PCM (\#152), Skew-Normal (\#98)], [Set $alpha >= eta_"out"$; PR-AUC on red-team injections; $o$-to-$p$ promotion lag],
    [*5. Quant. Finance & Customer Analytics*], [Mixed tabular, elliptical factors, warped financial time series], [One-hot $L_2$ distortion, singular $Sigma^(-1)$ inversion, lead-lag $L_2$ error], [K-Prototypes (\#8), FasterPAM (\#24), MCLUST (\#93), Ward (\#59), k-Shape (\#162), Soft-DTW (\#164)], [PAM (\#4), W-K-Means (\#28), MWK (\#29), DBA (\#163), ROCKET (\#165), Mixmod (\#94)], [BIC across 14 MCLUST shapes; Gower weight $gamma$; out-of-sample HRP drawdown],
    [*6. High-Velocity IoT & Edge Streaming*], [$n -> infinity$, edge RAM $< 1$ GB, anytime bursts, concept drift], [Multi-pass RAM overflow, stale history dominance, dropped packets], [BIRCH (\#64), BICO (\#78), StreamKM++ (\#195), ClusTree (\#196), DenStream (\#194), D-Stream (\#197)], [CluStream (\#193), Mini-Batch KM (\#3), BUBBLE (\#71), PERCH (\#76), GRINCH (\#77)], [Decay half-life $lambda$; coreset $epsilon$-budget; sliding-window SSE & drift alarms],
    [*7. Geospatial & Urban Mobility*], [Geodesic $(phi, lambda) times t times a$, urban/rural density ratios, GPS tracks], [Unit-scale failure of $(x, y, t)$ $L_2$, single-$epsilon$ merge, trajectory splits], [ST-DBSCAN (\#40), TRACLUS (\#169), OPTICS (\#33), VDBSCAN (\#50), WaveCluster (\#135)], [STING (\#133), CLARANS (\#6), GDBSCAN (\#39), GriDBSCAN (\#48), MR-DBSCAN (\#49)], [Decouple $(epsilon_1, epsilon_2, Delta epsilon)$; Haversine metric; MDL turning-point audit],
    [*8. Social, Web & Biological Networks*], [Sparse $G=(V, E)$, $|V| <= 10^9$, power-law degrees, overlapping roles], [Modularity resolution limit, SBM degree-grouping, disconnected Louvain], [Leiden (\#108), Infomap (\#113), DC-SBM (\#131), BigCLAM (\#125), OSLOM (\#124), METIS (\#121)], [Louvain (\#107), MCL (\#109), SBM (\#130), CPM (\#127), SLPA (\#126), DEMON (\#128)], [Check Kesten–Stigum SNR; CPM $gamma$ sweep; conductance, NMI, edge-cut balance],
    [*9. Industrial Sensors & Subspaces*], [$d in [10^2, 5 times 10^3]$, cluster-specific $d_j << d$, rotated sensor faults], [Ambient distance concentration, global PCA mixes disjoint fault axes], [SSC (\#150), ORCLUS (\#138), PROCLUS (\#136), STATPC (\#149), P3C (\#148), COSA (\#145)], [CLIQUE (\#134), SUBCLU (\#137), FIRES (\#146), INSCY (\#147), OptiGrid (\#141)], [Separate axis-parallel vs. oblique tests; $ell_1$ tolerance $lambda$; subspace F1],
    [*10. Multi-Modal Fusion & Consensus*], [$V >= 2$ heterogeneous views, mixed metrics, model-selection ambiguity], [Concatenation view-dominance, single-view noise, seed/family variance], [Co-reg. Spectral (\#175), Multi-View KM (\#174), MCLA (\#173), CSPA (\#171), VaDE (\#184), IDEC (\#182)], [HGPA (\#172), Monti (\#170), DEC (\#181), DCN (\#183), Kernel KM (\#16)], [ANMI consensus maximization; view weight audit; reconstruction loss (IDEC/VaDE)],
  )
]

#wide_figure(caption: [*Table II:* Observable Data-Invariant Regime to Provably Admissible Algorithm Families (Synthesizing Theorem 1, Rules 1--6). Identifies admissible candidates, disqualified families, operational complexity, and a recommended three-paradigm ClusterOverview verification trio.])[
  #set par(justify: false, first-line-indent: 0pt)
  #set text(size: 6.7pt)
  #show math.equation: set text(size: 6.7pt)
  #table(
    columns: (1.20fr, 1.38fr, 1.22fr, 0.85fr, 1.35fr),
    align: (left, left, left, left, left),
    inset: (x: 3.0pt, y: 2.5pt),
    fill: (_, y) => if y == 0 { rgb("#f1f5f9") } else if calc.odd(y) { rgb("#fbfdff") } else { white },
    table.header(
      [Observable Invariant Regime $cal(I)(cal(X))$], [Provably Admissible Families (Exact IDs)], [Eliminated Families (Theorem 1 Rule)], [Asymptotic Cost], [Recommended ClusterOverview Verification Trio],
    ),
    [*A. Clean Convex / Ellipsoidal* ($n <= 10^5, d <= 30, eta < 0.01$)], [K-Means++ (\#2), Hartigan–Wong (\#19), GMM (\#81), MCLUST (\#93), Ward (\#59), GK (\#153)], [Exhaustive subspace (Tier 6 overhead); streaming heuristics unnecessary], [$O(n k d T)$ to $O(n k d^2 T)$], [K-Means++ (\#2) + MCLUST/GMM (\#81/\#93) + Ward (\#59)],
    [*B. Massive Batch Euclidean* ($n > 10^6, d <= 200$, RAM-bound)], [Mini-Batch KM (\#3), K-Means|| (\#23), BIRCH (\#64), BICO (\#78), Elkan/Hamerly (\#13--14), CLARA (\#5)], [All Tier 4--5 $O(n^2)$ matrix methods (Rule 1: RAM overflow)], [$O(n d)$ or $O(b k d T)$], [Mini-Batch KM (\#3) + BIRCH (\#64) + BICO (\#78)],
    [*C. Unbounded Stream & Concept Drift* ($n = infinity, B_"mem" = o(n d)$)], [DenStream (\#194), ClusTree (\#196), D-Stream (\#197), CluStream (\#193), StreamKM++ (\#195)], [All batch multi-pass algorithms (Algs. 1--2, 4, 31--34, 56--63, 81--105) (Rule 1)], [$O(1)$--$O(q d)$ per event], [DenStream (\#194) + ClusTree (\#196) + StreamKM++ (\#195)],
    [*D. Non-Convex / Multi-Density Manifolds* ($d_"int" <= 15$)], [HDBSCAN (\#32), OPTICS (\#33), SNN (\#36), Spectral NJW/Ncut (\#103--104), Leiden-kNN (\#108), DPC (\#51)], [All convex Bregman partitioners (Cat. I); fixed-$epsilon$ DBSCAN if density varies (Rule 4)], [$O(n log n)$ to $O(n^2)$], [HDBSCAN (\#32) + Spectral NJW (\#103) + Leiden-kNN (\#108)],
    [*E. Heavy Outliers / Adversarial Noise* ($eta_"out" in [0.02, 0.30]$)], [HDBSCAN (\#32), OPTICS-OF (\#53), Trimmed KM (\#26), TMM (\#88), FasterPAM (\#24), K-Medians (\#15)], [Unweighted $1 / n$-breakdown means & Gaussian MLEs (Algs. 1, 59, 81, 151) (Rule 3)], [$O(n k d T)$ to $O(n^2)$], [TMM (\#88) + Trimmed KM (\#26) + HDBSCAN (\#32)],
    [*F. High-$d$ Sparse / Spherical Embeddings* ($d >= 300, rho_"sp" > 0.9$)], [Spherical KM (\#12), vMF (\#89), LDA/CTM (\#83/\#100), DMM/BTM (\#97/\#101), NMF (\#92), SNN (\#36)], [Ambient $L_2$ Voronoi, $k$-d tree DBSCAN, dense covariance GMM (Rules 2 & 5)], [$O("nnz"(X) k T)$], [Spherical KM (\#12) + vMF (\#89) + Leiden-cosine (\#108)],
    [*G. High-$d$ Embedded Subspaces* ($d > 100, d_j << d$)], [SSC (\#150), ORCLUS (\#138), PROCLUS (\#136), STATPC (\#149), P3C (\#148), COSA (\#145)], [Full-space $L_2$ (Rule 2); CLIQUE/SUBCLU if $d_j > 12$ (Ch. XII Thm. 1)], [$O(n k d)$ to $O(n^3)$], [PROCLUS (\#136) + ORCLUS (\#138) + SSC (\#150)],
    [*H. Mixed Numeric–Categorical Tabular* ($iota_"int" = "strict"$)], [K-Prototypes (\#8), FasterPAM/Gower (\#24), LCA (\#90), Mixmod (\#94), ROCK (\#67), AutoClass (\#87)], [Raw $L_2$ on one-hot codes (Rule 5); opaque deep encoders (Rule 6)], [$O(n k d T)$ to $O(n^2)$], [K-Prototypes (\#8) + FasterPAM (\#24) + Mixmod/LCA (\#90/\#94)],
    [*I. Sparse Relational / Overlapping Graphs* ($G = (V, E, W)$)], [Leiden (\#108), Infomap (\#113), DC-SBM (\#131), BigCLAM (\#125), OSLOM (\#124), METIS (\#121)], [Euclidean embedding without spectral/flow invariance (Rule 5); naive SBM on hubs], [$O(m)$ to $O(m log n)$], [Leiden (\#108) + Infomap (\#113) + DC-SBM (\#131)],
    [*J. Elastic Time Series & Multi-Modal Views*], [k-Shape (\#162), Soft-DTW (\#164), DBA (\#163), Multi-View KM (\#174), Co-reg. Spectral (\#175), MCLA (\#173)], [Pointwise $L_2$ (breaks phase/warp invariance); raw feature concatenation (Rule 5)], [$O(n k m log m)$ to $O(V n^2)$], [k-Shape (\#162) + Soft-DTW (\#164) + MCLA Consensus (\#173)],
  )
]

// =============================================================================
// SECTION IV: PRODUCTION ENGINEERING CHECKLIST & CLUSTEROVERVIEW VERIFICATION
// =============================================================================
#sec("IV", "Production Engineering Checklist, Validation Pitfalls, and ClusterOverview Verification")

Even an algorithm chosen rigorously from Tables I--II fails operationally if deployed atop flawed preprocessing, biased validation metrics, or unverified single runs. We conclude with the mathematical safeguards uniting theory with audited production engineering.

#subsec("A", "Preprocessing, Feature Scaling, Imputation, and Metric Hygiene")

Unlike supervised learners that absorb column scaling into regression coefficients, unsupervised clusterers take input distances literally. The next proposition establishes the exact threshold at which an unscaled noise feature hijacks 2-Means and isolates the affine-invariant model class.

#proposition(title: "Feature-Scaling Signal Hijacking and Affine-Invariance Dichotomy", origin: "Original Proof")[
  Let $X = (X_1, dots, X_d) in bb(R)^d$ combine two equal-weight clusters along coordinate $1$ ($X_1 = plus.minus Delta / 2 + Z_1$, $Z_1 tilde.op cal(N)(0, 1)$, so $"Var"(X_1) = 1 + Delta^2 / 4$) with independent pure Gaussian noise coordinates $X_r tilde.op cal(N)(0, kappa)$ ($r = 2, dots, d$, $kappa > 1$).
  1. *(Isotropic $L_2$ Signal Hijacking):* Under unscaled Euclidean distance (Algs. 1--3, 15, 31, 56--59, 151), signal-to-noise degrades by $O(1 / kappa)$:
     $ "SNR"_"raw" = Delta^2 / (2 (1 + (d - 1) kappa)) = O( "SNR"_"std" / kappa ), $
     and bisecting a pure-noise coordinate $X_r = 0$ achieves strictly lower 2-Means SSE (Alg. 1) than *any* partition along the true signal coordinate $X_1$ whenever
     $ kappa >= (pi / 2) (1 + Delta^2 / 4) approx 1.571 (1 + Delta^2 / 4). $
  2. *(Exact Affine Invariance):* Full-covariance GMM / MCLUST-VVV (Algs. 81, 93 @fraley2002mclust), full-$Sigma$ TMM (Alg. 88 @peel2000robust), and Gustafson–Kessel (Alg. 153 @gustafson1979gk) are *strictly invariant* under every invertible affine map $x |-> A x + b$ ($A in "GL"(d, bb(R))$).
]

#proof[
  *(1)* Any 2-partition determined purely by coordinate $X_1$ leaves all noise variances untouched and can reduce coordinate-1 variance from $"Var"(X_1) = 1 + Delta^2 / 4$ by *strictly less* than $1 + Delta^2 / 4$ (since within-cluster Gaussian variance remains strictly positive). Conversely, bisecting pure-noise coordinate $X_2 tilde.op cal(N)(0, kappa)$ at $0$ yields symmetric half-normal centroids $plus.minus bb(E)[|X_2|] = plus.minus sqrt(2 kappa / pi)$, achieving exact SSE reduction $(bb(E)[|X_2|])^2 = (2 / pi) kappa$. Therefore, splitting the pure noise axis strictly beats every signal-axis partition whenever $(2 / pi) kappa >= 1 + Delta^2 / 4$, i.e., $kappa >= (pi / 2)(1 + Delta^2 / 4)$.
  *(2)* Under $tilde(x)_i = A x_i + b$ ($det A != 0$) with $(tilde(mu)_j, tilde(Sigma)_j) = (A mu_j + b, A Sigma_j A^top)$,
  $ & (tilde(x)_i - tilde(mu)_j)^top tilde(Sigma)_j^(-1) (tilde(x)_i - tilde(mu)_j) \
    &quad = (x_i - mu_j)^top A^top (A^(-top) Sigma_j^(-1) A^(-1)) A (x_i - mu_j) \
    &quad = (x_i - mu_j)^top Sigma_j^(-1) (x_i - mu_j). $
  The constant Jacobian $|det A|^(-1)$ cancels from all E-step responsibilities $gamma_(i j)$ in GMM/MCLUST/TMM, and $|det A|^(2 / d)$ cancels identically from the determinant-normalized fuzzy covariance norm in Gustafson–Kessel (\#153).
]

*Mandatory Preprocessing Hygiene Rules:*
1. *ID & Constant Pruning:* Drop surrogate keys, integer timestamps, and zero-variance columns (automated in ClusterOverview's `src/data.rs`) to prevent artificial 1D line manifolds.
2. *Distribution-Matched Scaling:* Use Z-score only on light-tailed features; use *IQR/MAD robust scaling* when $eta_"out" > 0.01$ (preventing outliers from crushing clean variance); apply *log1p / CLR* to counts/compositions; and *row-wise $L_2$ normalization* on embeddings.
3. *Missing-Data Hygiene:* Avoid mean/zero imputation at $> 15\%$ missingness (which injects Dirac density spikes creating false DBSCAN/GMM clusters); prefer iterative/kNN imputation (MAR) or missingness flags under Gower/K-Prototypes (MNAR).
4. *PCA Denoising Without Blind Whitening:* Truncated PCA ($80\%$--$95\%$ variance) suppresses ambient noise before SNN/Leiden or HDBSCAN, whereas full whitening ($Sigma^(-1 / 2)$) inflates trailing noise eigenvectors to unit scale, reinstating Rule 2 collapse.

#subsec("B", "Seven Recurring Production Validation Pitfalls and Mathematical Remedies")

1. *Pitfall 1 — Convexity & Noise-Trimming Bias of Internal Indices:* Silhouette @rousseeuw1987silhouettes, CHI @calinski1974dendrite, DBI @davies1979cluster, and Dunn @dunn1974well systematically penalize non-convex manifolds and reward aggressive noise excision (Chapter XI, Section VI). In ClusterOverview's Iris benchmark (Chapter XIII), *OPTICS* (\#33) and *MAFIA* (\#140) achieve top Silhouette ($0.657, 0.642$) by discarding $21\%$--$25\%$ of boundary samples as noise and collapsing 3 species into $k = 2$ ($"ARI" = 0.442, 0.388$), whereas *GMM* (\#81) keeps 100% coverage at $k = 3$ and achieves $"ARI" = 0.950$. _Remedy:_ Always co-report noise ratio, balance entropy, density validity (DBCV), or graph conductance alongside Silhouette.
2. *Pitfall 2 — Clustering on 2D UMAP / t-SNE Embeddings:* Nonlinear 2D projections intentionally equalize local density and distort global distances—invalidating DBSCAN/HDBSCAN and fabricating gaps in uniform noise. _Remedy:_ Cluster in native, truncated-PCA ($d_"pca" in [15, 50]$), or SNN graph space; reserve 2D embeddings strictly for visual inspection.
3. *Pitfall 3 — The $k$-Selection Mirage & Stability Paradox:* Ben-David et al. @ben-david2006sober proved that symmetric true $k$-cluster geometries can exhibit *lower* bootstrap stability than wrong $k' != k$. _Remedy:_ Triangulate $k$ via penalized likelihood (BIC/ICL in MCLUST \#93), topological persistence (HDBSCAN \#32 / Monti CDF \#170 @monti2003consensus), and downstream domain utility.
4. *Pitfall 4 — Silent $O(n^2)$ Spatial-Tree Degradation:* Above $d_"int" > 15$, $k$-d and ball trees lose pruning power, silently turning DBSCAN/OPTICS into $O(n^2)$ scans @schubert2017dbscanrevisited. _Remedy:_ Monitor tree pruning ratios and transition to HNSW/IVF $k$-NN graphs + Leiden (\#108) or Cat. VI subspaces.
5. *Pitfall 5 — Seed Sensitivity & Row-Order Dependence:* Across $cal(A)_(200)$, 101 algorithms are seed/stochastic and 31 are row-order dependent (BIRCH \#64, ART \#178, COBWEB \#85). _Remedy:_ Randomly shuffle sorted SQL extracts before tree ingestion, use K-Means++ (\#2 @arthur2007kmeanspp) with $R >= 20$ seeded restarts, and pin RNG seeds.
6. *Pitfall 6 — Resolution-Limit Over-Merging & Deep Collapse:* Default modularity merges small cliques below $sqrt(m / 2)$ edges (Theorem 4 @fortunato2007resolution); unregularized deep clusterers collapse onto 1--2 prototypes. _Remedy:_ Audit `balance_entropy`, sweep CPM $gamma$, and enforce Sinkhorn/reconstruction losses (IDEC \#182, VaDE \#184, SwAV \#190).
7. *Pitfall 7 — Preprocessing Leakage & Unmonitored Drift:* Fitting scalers across train/test splits leaks distribution shift; static models silently degrade under live covariate drift. _Remedy:_ Freeze versioned preprocessing pipelines and track live Mahalanobis drift and `noise_ratio` spikes in production.

#subsec("C", "Multi-Algorithm Ensemble Verification via ClusterOverview")

Without ground-truth labels, engineering confidence requires *Cross-Paradigm Agreement Verification*: comparing independent partitions from orthogonal algorithm families via chance-corrected Adjusted Rand Index @hubert1985comparing.

#theorem(title: "Cross-Paradigm Agreement Certificate Against Inductive-Bias Artifacts", origin: "Original Proof")[
  Let $Pi_1 = A_1 (cal(X))$ and $Pi_2 = A_2 (cal(X))$ be $k$-cluster partitions from two structurally distinct families ($n_(1, a), n_(2, b) <= c_"max" n < n$).
  1. *(Null Concentration Under $cal(H)_0$):* Under random hypergeometric alignment between $Pi_1$ and $Pi_2$ (structureless null $cal(H)_0$), $bb(E)_0 ["ARI"(Pi_1, Pi_2)] = 0$,
     $ "Var"_0 ("ARI"(Pi_1, Pi_2)) = O( 1 / n ), \
       bb(P)_0 ("ARI"(Pi_1, Pi_2) >= tau) <= O( 1 / (n tau^2) ). $
  2. *(Triangle Lower Bound to Unknown Ground Truth):* Letting $d_"RD" (Pi, Pi') := 1 - "RI"(Pi, Pi')$ denote pair-disagreement distance, for *every* unknown true partition $Pi^*$,
     $ max( d_"RD" (Pi_1, Pi^*), d_"RD" (Pi_2, Pi^*) ) >= (1 / 2) d_"RD" (Pi_1, Pi_2). $
     Thus low cross-paradigm ARI provably certifies severe misspecification in at least one family, whereas high cross-paradigm ARI ($>= 0.75$) across orthogonal families rejects $cal(H)_0$ at $O(1 / n)$ significance.
]

#proof[
  *(1)* With contingency counts $N_(a b) = |C_(1, a) inter C_(2, b)|$, $T = sum_(a, b) binom(N_(a b), 2)$, $P = sum_a binom(n_(1, a), 2)$, $Q = sum_b binom(n_(2, b), 2)$, and $M = binom(n, 2)$,
  $ "ARI"(Pi_1, Pi_2) = (T - (P Q) / M) / ((1 / 2)(P + Q) - (P Q) / M) $
  has exact expectation $0$ under hypergeometric alignment @hubert1985comparing. Decomposing $"Var"_0 (T)$ over identical, vertex-sharing, and disjoint edge pairs yields $"Var"_0 (T) = O(n^3)$ against squared denominator $Theta(n^4)$, giving $"Var"_0 ("ARI") = O(1 / n)$ and Chebyshev bound $O(1 / (n tau^2))$.
  *(2)* Embedding partitions as binary pair vectors $v(Pi) in {0, 1}^M$ makes $d_"RD" (Pi, Pi') = M^(-1) norm(v(Pi) - v(Pi'))_1$ an exact metric; the triangle inequality gives
  $ d_"RD" (Pi_1, Pi_2) &<= d_"RD" (Pi_1, Pi^*) + d_"RD" (Pi^*, Pi_2) \
    &<= 2 max( d_"RD" (Pi_1, Pi^*), d_"RD" (Pi_2, Pi^*) ), $
  proving the lower bound immediately.
]

*Operationalizing Verification in ClusterOverview:* The pure-Rust *ClusterOverview* desktop and CLI platform operationalizes this workflow natively across all 200 algorithms: (1) *Shared Precomputation (`src/data.rs`):* decodes CSVs, prunes ID/constant columns, computes Z-score/Min-Max scaling, 2D PCA, pairwise distances, and $k$-NN graphs once into read-only `Arc<PreparedData>`; (2) *Parallel Multi-Paradigm Sweep (`src/algorithms/`):* executes shortlists or the full 200-algorithm catalog across CPU cores via `rayon` with isolated error recovery; (3) *Seven-Metric Audit (`src/metrics.rs`):* scores Silhouette, DBI, CHI, Dunn, SSE, Balance Entropy, and Noise Ratio simultaneously, exposing Pitfalls 1 and 6 immediately; (4) *Interactive A-vs-B ARI Inspection (`src/gui.rs`):* verifies Table II cross-paradigm trios quantitatively via Theorem 6 alongside 2D scatter overlays; and (5) *Consensus & Audit Export (`src/pdf_report.rs`):* promotes the Pareto-optimal survivor (or MCLA/CSPA consensus \#173/\#171 @strehl2002cluster) and generates a reproducible vector-PDF governance dossier (GUI or headless CLI `--csv ... --pdf ...`). Table III codifies the complete six-stage sign-off checklist.

#wide_figure(caption: [*Table III:* End-to-End Production Clustering Engineering and Verification Checklist. Unites formal theory (Theorems 1--6, Propositions 1--8) with operational verification in ClusterOverview and live production monitoring.])[
  #set par(justify: false, first-line-indent: 0pt)
  #set text(size: 6.7pt)
  #show math.equation: set text(size: 6.7pt)
  #table(
    columns: (0.95fr, 1.85fr, 1.55fr, 1.55fr),
    align: (left, left, left, left),
    inset: (x: 3.0pt, y: 2.5pt),
    fill: (_, y) => if y == 0 { rgb("#f1f5f9") } else if calc.odd(y) { rgb("#fbfdff") } else { white },
    table.header(
      [Lifecycle Stage], [Mandatory Technical Checks & Invariants], [Failure Mode Prevented], [Verification Artifact / Metric],
    ),
    [*1. Ingestion & Profiling*], [Drop surrogate IDs, raw timestamps, zero-variance cols; audit missingness (MAR/MNAR), sparsity $rho_"sp"$, tails, mixed types; estimate $(n, d, d_"int", eta_"out")$], [Spurious ID line manifolds; silent type coercion; missing-value Dirac density spikes], [Column profile (`data.rs`); kurtosis & outlier quantiles; signature $cal(I)(cal(X))$],
    [*2. Scaling & Metric Design*], [Match transform to metric invariance: Z-score (light-tail), IQR/MAD (heavy-tail), log1p/CLR (counts), $L_2$-unit ($bb(S)^(d-1)$), Gower/Haversine/SBD; forbid blind PCA whitening], [Scale hijacking (Prop. 8); unit collapse on $(x, y, t)$ (Prop. 6); length confounding], [Post-scaling variance & condition audit; serialized versioned scaler],
    [*3. Selection Funnel (Alg. 1)*], [Apply Theorem 1 Rules 1--6 via Stages 1--4: enforce RAM/streaming, high-$d$, breakdown, non-convexity, and $iota_"int"$ gates; shortlist $3$--$6$ cross-paradigm candidates], [OOM crashes on $O(n^2)$ matrices; $1 / n$ outlier breakdown; convex bias; opaque models], [Elimination log citing Rules 1--6; Table II shortlist $cal(S)_4$],
    [*4. Execution & Stability*], [Shuffle row order for order-dependent trees (BIRCH/ART); pin seeds; run $R >= 20$ K-Means++ restarts; audit $k$-d pruning for $d_"int" > 15$; tune CPM $gamma$ or $alpha >= eta_"out"$], [Row-order tree skew; poor local minima; silent $O(n^2)$ tree collapse; resolution over-merging], [Multi-seed objective variance; Monti CDF stability curve; latency/RAM profile],
    [*5. ClusterOverview Audit*], [Run shortlist (or 200-alg sweep) via `rayon` on `Arc<PreparedData>`; audit all 7 metrics simultaneously; compute pairwise A-vs-B cross-paradigm ARI; test MCLA/CSPA], [Silhouette noise-trimming illusion (Pitfall 1); single-family bias (Thm. 6); entropy collapse], [ClusterOverview vector-PDF dossier; cross-paradigm ARI matrix ($>= 0.70$); noise/entropy gate],
    [*6. Production & Drift Ops*], [Freeze versioned scalers & prototypes; decouple online inference from offline refits; monitor live Mahalanobis/Wasserstein drift, `noise_ratio`, and balance entropy], [Training-serving skew; silent concept drift; gradual cluster extinction or outlier floods], [Live drift alerts; automated re-clustering SLA triggers; downstream KPI audit],
  )
]

// =============================================================================
// CONCLUDING MONOGRAPH SYNTHESIS
// =============================================================================
#sec("V", "Concluding Synthesis: Unifying Mathematical Foundations and Engineering Practice")

Across Parts I through IX of this monograph, two hundred clustering algorithms spanning seven decades of statistical, geometric, spectral, topological, and neural research have been unified onto a single formal mathematical structure: six variational-geometric lenses (Chapter I), exact asymptotic reductions (Chapter II), sharp axiomatic impossibility and concentration ceilings (Chapter XI), six stratified complexity tiers (Chapter XII), empirical multi-geometry diagnostics (Chapter XIII), and the ten-vertical engineering architecture completed in this chapter.

The decisive synthesis for theorists and practitioners alike is unmistakable: *rigorous unsupervised learning is never the pursuit of a single universal clustering algorithm, but the disciplined alignment of algorithmic invariants with the geometry, noise physics, computational constraints, and epistemic goals of the data.* By enforcing Theorem 1's elimination calculus, honoring the domain and preprocessing invariants of Propositions 2--8, avoiding the seven validation pitfalls of Section IV-B, and certifying partitions via cross-paradigm ARI agreement in ClusterOverview (Theorem 6), practitioners elevate clustering from fragile exploratory heuristics into a provably grounded, reproducible engineering discipline.
