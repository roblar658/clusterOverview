#import "../template/ieee.typ": *

#cmp_entry(
  roman: "XII",
  title: "Complete 200-Algorithm Master Comparative Reference Matrices",
  scope_label: "Exhaustive 10-Dimension Theoretical Reference Compendium for Algorithms 1-200",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)

#sec("I", "Architecture of the Ten-Dimension Reference Compendium")

Comparison Chapters I–XI have developed the grand variational unification across all eight clustering families, proved fifteen cross-family algebraic equivalences, dissected intra-category trade-offs, and established the axiomatic impossibility and robustness frontiers governing unsupervised learning @jain2010data, @xu2005survey, @aggarwal2013data. To transform those theoretical insights into an actionable desk reference for applied mathematicians, statisticians, and systems engineers, this chapter compiles the *complete, zero-gap, 200-row Master Comparative Reference Compendium* covering every single algorithm from Alg. 1 through Alg. 200 across *ten standardized theoretical dimensions*.

Rather than presenting a flat catalog without analytical structure, we organize Chapter XII into two tightly coupled halves. Sections I–IV establish the *macro-structural laws* that emerge when all 200 algorithms are projected onto complexity, determinism, memory-footprint, and high-dimensional axes—proving three formal theorems (Theorems 1–3), one sensitivity proposition (Proposition 1), and a ten-axis Pareto conservation corollary (Corollary 1), supported by three dynamically reconciled cross-category synthesis tables (Tables I–III). Section V then switches cleanly to a full-width IEEE reference layout (`columns: 1`) to present all eight exhaustive booktabs master matrices (Tables IV–XI, 200 rows total), before Section VI restores standard two-column IEEE formatting to synthesize cross-category Pareto archetypes and bridge directly to the 800-run empirical benchmarks in Comparison Chapter XIII.

#subsec("A", "Standardized Ten-Dimension Schema and Asymptotic Notation")

Every row in Tables IV–XI records ten exact, algorithm-specific attributes verified directly against the formal derivations in Parts I–VIII:
1. *ID:* Canonical catalog index $1, dots, 200$ matching the chapter numbering of Parts I–VIII.
2. *Algorithm (Year):* Authoritative monograph name and landmark publication year.
3. *Time Complexity:* Tight operational worst-case or indexed bound in unified notation: $n$ (or $N$) observations, $d$ ambient dimension, $k$ (or $K$) clusters, $T$ (or $I, E$) iterations/epochs, $m$ graph edges (Cat. V) or time-series length (Cat. VII), $b$ mini-batch size, $s$ sample/coreset size, $V$ vocabulary size, and $"nnz"(X)$ non-zeros.
4. *Space:* Peak working-memory footprint excluding read-only external disk streams.
5. *Det. (Determinism Class):* Strictly four mutually exclusive regimes: `Det.` (strictly deterministic and permutation-invariant), `Order` (deterministic given a fixed stream or tie-breaking order, but permutation-sensitive), `Seed` (deterministic given a random initialization or subsample), and `Stoch.` (intrinsically stochastic dynamics such as mini-batch SGD, Gibbs MCMC, randomized asynchronous label sweeps, or evolutionary metaheuristics).
6. *Geometry / Bias:* Native decision-boundary geometry, metric manifold, statistical family, or structural invariant.
7. *Noise / Outlier:* Explicit outlier excision ($C_(-1)$ noise set, $alpha$-trimming, homeless nodes), breakdown-point class ($1 / n$ vs. medoid / $L_1$ / bounded-influence heavy tails), or chaining/collapse vulnerability.
8. *Auto-$k$:* Whether cluster cardinality $k$ is inferred automatically (`Yes` with mechanism: BIC, AD-test, $epsilon$-reachability, EoM, DP/CRP, modularity $Q$, MDL, vigilance $rho$, grid cuts), produced across all scales (`Dendro.`), or fixed externally (`No`).
9. *High-$d$:* Concentration resilience as $d -> infinity$ @beyer1999nearest, stratified into `Poor`, `Fair`, `Good`, and `Strong` with the surviving mathematical mechanism noted concisely.
10. *Primary Domain:* Signature scientific, biomedical, geospatial, network, or industrial workload where the algorithm's inductive bias aligns with real data geometry.

#wide_figure(caption: [*Table I:* Cross-Category Executive 10-Dimension Profile Across Categories I–VIII (Algorithms 1–200). Counts reconcile strictly to the 200-algorithm master matrices in Tables IV–XI; _Auto-$k$_ counts emergent, model-selected, or multi-resolution dendrogram methods; _High-$d$ (G/S)_ counts algorithms rated _Good_ or _Strong_ under high-dimensional concentration.])[
  #set par(justify: false, first-line-indent: 0pt)
  #set text(size: 7.0pt)
  #show math.equation: set text(size: 7.0pt)
  #table(
    columns: (1.55fr, 0.36fr, 0.78fr, 0.78fr, 0.92fr, 0.54fr, 0.60fr, 0.64fr, 0.96fr),
    align: (left, center, left, left, center, center, center, center, left),
    inset: (x: 3.2pt, y: 2.8pt),
    fill: (_, y) => if y == 0 { rgb("#f1f5f9") } else if y == 9 { rgb("#eff6ff") } else if calc.odd(y) { rgb("#fbfdff") } else { white },
    table.header(
      [Category (IDs)], [$N$], [Modal Time Tier], [Modal Space Tier], [Det / Ord / Seed / Sto], [Auto-$k$], [Noise Rob.], [High-$d$ G/S], [Dominant Inductive Bias],
    ),
    [*Cat. I:* Centroid (1–30)], [30], [$O(n k d T)$ (21)], [Linear $O(n d)$ (21)], [4 / 2 / 22 / 2], [3/30], [10/30], [6/30], [Convex Bregman / metric Voronoi],
    [*Cat. II:* Density & Mode-Seeking (31–55)], [25], [$O(n log n)$ (12)], [Linear $O(n d)$ (17)], [16 / 7 / 2 / 0], [25/25], [24/25], [3/25], [Non-convex KDE / $epsilon$-level sets],
    [*Cat. III:* Hierarchical & Tree-Based (56–80)], [25], [$O(n^2)$ (13)], [Quad. $O(n^2)$ (13)], [19 / 5 / 1 / 0], [19/25], [14/25], [7/25], [Ultrametric trees & CF coresets],
    [*Cat. IV:* Distribution (81–102)], [22], [$O(n k d T)$ (20)], [Bilinear $O(n k)$ (11)], [1 / 2 / 14 / 5], [10/22], [6/22], [13/22], [Exp.-family / topic / Bayes mix.],
    [*Cat. V:* Graph (103–132)], [30], [$O(n)$ / stream (8)], [Linear $O(n d)$ (18)], [12 / 4 / 5 / 9], [19/30], [13/30], [29/30], [Laplacian cuts, flow, SBM blocks],
    [*Cat. VI:* Grid (133–150)], [18], [Exp. / output (5)], [Linear $O(n d)$ (8)], [11 / 1 / 4 / 2], [14/18], [18/18], [15/18], [Axis / oblique low-$d$ subspaces],
    [*Cat. VII:* Fuzzy (151–175)], [25], [$O(n k d T)$ (13)], [Bilinear $O(n k)$ (15)], [5 / 3 / 15 / 2], [7/25], [10/25], [14/25], [Fuzzy simplex, 2-mode, TS warp],
    [*Cat. VIII:* Neural (176–200)], [25], [$O(n k d T)$ (15)], [Cov. / Deep (12)], [0 / 7 / 2 / 16], [5/25], [8/25], [15/25], [Neural latent, stream, metaheur.],
    [*Total Across Monograph (1–200)*], [*200*], [*Bilinear (80, 40%)*], [*Linear (78, 39%)*], [*68 / 31 / 65 / 36*], [*102 (51%)*], [*103 (51.5%)*], [*102 (51.0%)*], [*6 Master Lenses (Ch. I)*],
  )
]

Inspecting Table I reveals four striking global invariants of the 200-algorithm literature:
- *Exact Macro-Symmetry in Model-Order Selection:* Almost exactly half of all clustering methods require an external cluster count (98 / 200, 49.0%), while 102 / 200 (51.0%) infer $k$ automatically or output a multi-scale dendrogram—dominated by Category II (25 / 25, 100%), Category III (19 / 25, 76.0%), Category VI (14 / 18, 77.8%), and Category V (19 / 30, 63.3%).
- *Near-Equal Determinism Split:* Across all 200 methods, 99 (49.5%) are deterministic or order-deterministic (`Det.` = 68, `Order` = 31), while 101 (50.5%) rely on stochastic seeding or randomized optimization (`Seed` = 65, `Stoch.` = 36).
- *High-Dimensional Bifurcation:* Exactly 102 / 200 (51.0%) methods achieve _Good_ or _Strong_ high-dimensional resilience—concentrated heavily in sparse graphs (Cat. V: 29 / 30), subspace projectors (Cat. VI: 15 / 18), topic/directional models (Cat. IV: 13 / 22), and deep neural encoders (Cat. VIII: 15 / 25)—whereas 98 / 200 (49.0%) degrade under isotropic Euclidean distance concentration (notably Cat. I: 24 / 30 and Cat. II: 22 / 25).
- *Explicit Noise Protection:* 103 / 200 (51.5%) algorithms incorporate explicit outlier excision, heavy-tailed attenuation, medoid/median protection, or noise-robust regularization, while the remaining 97 methods retain the classical $1 / n$ arithmetic-mean breakdown or single-bridge chaining vulnerability.

#sec("II", "Asymptotic Complexity-Class Stratification Across 200 Algorithms")

Although textbook surveys routinely divide clustering algorithms into "fast" ($k$-means) and "slow" (hierarchical/spectral) methods, stratifying all 200 algorithms in Parts I–VIII exposes *six sharply separated asymptotic complexity tiers* spanning linear single-pass streaming ($O(n)$) to exponential subspace-lattice enumeration ($O(2^(k_"max") n d)$). Table II summarizes the exact population counts, typical space footprints, determinism rates, and canonical algorithm IDs across all six tiers.

#wide_figure(caption: [*Table II:* Six-Tier Asymptotic Complexity-Class Stratification Across Algorithms 1–200. Every algorithm is assigned to its canonical operational tier matching Parts I–VIII (exact counts sum strictly to 200). _Det.+Ord._ reports the deterministic/order-deterministic share within each complexity class.])[
  #set par(justify: false, first-line-indent: 0pt)
  #set text(size: 7.0pt)
  #show math.equation: set text(size: 7.0pt)
  #table(
    columns: (1.35fr, 1.15fr, 0.44fr, 0.90fr, 0.56fr, 1.35fr, 1.30fr),
    align: (left, left, center, left, center, left, left),
    inset: (x: 3.2pt, y: 2.8pt),
    fill: (_, y) => if y == 0 { rgb("#f1f5f9") } else if y == 7 { rgb("#eff6ff") } else if calc.odd(y) { rgb("#fbfdff") } else { white },
    table.header(
      [Complexity Tier], [Canonical Time Scaling], [Count (%)], [Typical Space], [Det.+Ord.], [Structural Bottleneck / Regime], [Representative Algorithm IDs],
    ),
    [*Tier 1: Linear & Sublinear Streaming*], [$O(n)$, $O(n d)$, $O(b k d T)$, $O(L m)$], [32 (16.0%)], [$O(q d)$ to $O(n + m)$], [19/32 (59%)], [Single-pass streams, sparse CSR graphs, subsampled medoids, fixed grids], [Algs. 3, 5, 23, 25, 47, 55, 64, 71, ..., 194, 195, 197],
    [*Tier 2: Log-Linear & Tree-Indexed*], [$O(n log n)$, $O(n d log n)$, $O(m log n)$], [28 (14.0%)], [$O(n d)$ or $O(n + m)$], [24/28 (86%)], [Spatial-index density, ANN/tree hierarchies, multilevel graph cuts], [Algs. 9, 31, 32, 33, 35, 37, 38, 39, ..., 167, 168, 196],
    [*Tier 3: Bilinear Iterative & Deep SGD*], [$O(n k d T)$, $O(n k d^2 T)$, $O(E n C_"net")$], [80 (40.0%)], [$O(n (d + k))$ or $O(|theta| + B d)$], [6/80 (8%)], [Prototype $k$-means, EM mixtures, topic models, fuzzy, deep SGD, bio-evo.], [Algs. 1, 2, 6, 7, 8, 10, 11, 12, ..., 198, 199, 200],
    [*Tier 4: Quadratic Pairwise & Kernel*], [$O(n^2)$, $O(n^2 d)$, $O(n^2 log n)$], [41 (20.5%)], [$O(n^2)$ (or $O(n d)$ stream)], [34/41 (83%)], [NN-chain linkages, PAM/FasterPAM, Kernel K-Means/FCM, Mean Shift, DPC], [Algs. 4, 16, 17, 24, 34, 36, 41, 43, ..., 170, 171, 186],
    [*Tier 5: Cubic Spectral & Exact Tree*], [$O(n^3)$, $O(m^2 n)$, $O(n^2 (d + n))$], [12 (6.0%)], [$O(n^2)$ or $O(V n^2)$], [10/12 (83%)], [Dense spectral eigensolvers, DIANA, Minimax, NJ, Girvan–Newman, SSC], [Algs. 42, 63, 73, 75, 103, 104, 105, 112, 118, 123, 150, 175],
    [*Tier 6: Output-Sensitive / Exponential*], [$O(2^(k_"max") n d)$, $O(d n 3^(d / 3))$], [7 (3.5%)], [$O(n d + |cal(U)|)$], [6/7 (86%)], [Exhaustive Apriori subspace lattices, exact $k$-cliques, binary Bimax], [Algs. 127, 134, 137, 139, 140, 143, 159],
    [*Total Monograph Compendium*], [—], [*200 (100%)*], [—], [*99/200 (49.5%)*], [*Full spectrum $O(n)$ to $O(2^s)$*], [*Complete IDs 1–200 (zero gaps)*],
  )
]

#subsec("A", "Quantitative Anatomy of the Six Complexity Tiers")

1. *Tier 1 — Linear & Sublinear Streaming (32 algs., 16.0%):* Anchored by single-pass additive-summary engines (BIRCH Alg. 64 @zhang1996birch, BICO Alg. 78 @fichtenberger2013bico, CluStream Alg. 193 @aggarwal2003clustream, DenStream Alg. 194 @cao2006denstream, #box[StreamKM++] Alg. 195 @ackermann2012streamkmpp, D-Stream Alg. 197 @chen2007dstream), subsampled/mini-batch quantizers (Mini-Batch K-Means Alg. 3 @sculley2010webscale, CLARA Algs. 5 & 25, #box[DBSCAN++] Alg. 47 @jang2019dbscanpp), fixed-cell spatial grids (STING Alg. 133 @wang1997sting, WaveCluster Alg. 135 @sheikholeslami1998wavecluster, SLIC Alg. 166 @achanta2012slic), and sparse-graph local movers (Louvain Alg. 107 @blondel2008louvain, Leiden Alg. 108 @traag2019leiden, LPA Alg. 110 @raghavan2007lpa, PIC Alg. 132 @lin2010pic). Every Tier 1 method avoids $Omega(n^2)$ pairwise materialization by restricting interaction to either $O(k)$ global summaries, $O(1)$ local grid stencils, or $O("deg"(v))$ sparse graph neighbors.
2. *Tier 2 — Log-Linear & Tree-Indexed (28 algs., 14.0%):* Dominated by spatial-tree density methods (DBSCAN, HDBSCAN, OPTICS, DENCLUE, DeBaCl, DBSCAN-3D, Algs. 31–33, 35, 46, 168 @ester1996dbscan, @campello2013hdbscan, @ankerst1999optics), incremental/ANN hierarchies (PERCH, GRINCH, SCC, FINCH, Algs. 76–77, 79–80), and priority-queue/multilevel graph engines (Infomap Alg. 113, CNM Alg. 115, METIS Alg. 121, Chaco Alg. 122). Crucially, *spatial-index Tier 2 bounds are dimension-conditional*: when intrinsic dimension $d_"int" > 15$, $k$-d trees and ball trees lose pruning efficiency, silently promoting Algs. 31–33, 38, 40, and 53 into Tier 4 ($O(n^2)$) @schubert2017dbscanrevisited.
3. *Tier 3 — Bilinear Iterative Prototype, Mixture & Deep SGD (80 algs., 40.0%):* The structural centerpiece of the entire monograph, accounting for exactly two-fifths of all 200 methods. It unites alternating Bregman partitioners ($O(n k d T)$ in Cat. I), parametric EM and variational topic models ($O(T n k d)$ to $O(T n k d^2 + T k d^3)$ in Cat. IV @dempster1977em, @blei2003lda), fuzzy and time-series centroid engines (Cat. VII @bezdek1981fcm, @paparrizos2015kshape, @cuturi2017softdtw), deep mini-batch representation learners ($O(E n (C_"net" + k p))$ in Cat. VIII @xie2016dec, @caron2020swav), and population metaheuristics ($O(G P n k d)$, Algs. 198–200). Notice that Tier 3 has the *lowest determinism share* (just 6 / 80, or 8% deterministic/order-deterministic), because bilinear alternating minimization over non-convex landscapes universally relies on randomized initialization or stochastic mini-batch gradients.
4. *Tier 4 — Quadratic Pairwise, Kernel & Agglomerative (41 algs., 20.5%):* Encompasses exact Lance–Williams NN-chain linkages (Algs. 56–62, 69–70, 72, 74 @mullner2011modern), exact metric medoid search (PAM Alg. 4, FasterPAM Alg. 24 @schubert2019fasterpam), RKHS kernel clusterers (Kernel K-Means Alg. 16 @dhillon2004kernel, KFCM Alg. 155 @zhang2004kfcm), exact unindexed mode seekers (Mean Shift Alg. 34, SNN Alg. 36, Quickshift Alg. 41, DPC Alg. 51 @rodriguez2014dpc), and co-association ensembles (Monti Alg. 170, CSPA Alg. 171 @strehl2002cluster). Conversely to Tier 3, Tier 4 exhibits a *high determinism share* (34 / 41, or 83% deterministic/order-deterministic), reflecting exact all-pairs graph or matrix inspection.
5. *Tier 5 — Cubic Spectral, Divisive & Exact Combinatorial (12 algs., 6.0%):* Contains twelve exact eigensolver, divisive, or optimization methods: Medoidshift (Alg. 42), DIANA (Alg. 63), Minimax Linkage (Alg. 73 @bien2011minimax), Neighbor-Joining (Alg. 75 @saitou1987neighbor), dense spectral cuts (NJW Alg. 103, Ncut Alg. 104, Ratio Cut Alg. 105 @vonluxburg2007tutorial), Girvan–Newman (Alg. 112), HCS (Alg. 118), Kluger Biclustering (Alg. 123), Sparse Subspace Clustering (Alg. 150 @elhamifar2013ssc), and Co-regularized Multi-View Spectral (Alg. 175 @kumar2011coregularized). Sparse Lanczos/ARPACK solvers reduce spectral iterations from $O(n^3)$ to $O(I_"Lan" m k)$ when $W$ is a sparse $k$-NN graph, yet dense affinity matrices remain strictly cubic.
6. *Tier 6 — Output-Sensitive & Exponential Combinatorial (7 algs., 3.5%):* Contains seven exhaustive combinatorial searchers: Clique Percolation Method (Alg. 127 @palla2005cpm), bottom-up Apriori subspace finders (CLIQUE Alg. 134, SUBCLU Alg. 137, DOC Alg. 139, MAFIA Alg. 140, ENCLUS Alg. 143), and exact binary biclustering Bimax (Alg. 159 @prelic2006bimax). Every Tier 6 algorithm trades polynomial worst-case guarantees for *combinatorial completeness* over overlapping cliques, subspaces, or inclusion-maximal submatrices.

#subsec("B", "Fundamental Information-Theoretic Complexity Barriers")

Why do these six complexity tiers persist after sixty years of algorithmic engineering? The next theorem proves that the sharp jumps between Tier 1 ($O(n d)$), Tier 4 ($O(n^2)$), and Tier 6 ($O(2^s)$) reflect irreducible information-theoretic barriers rather than suboptimal implementations.

#theorem(title: "Information-Theoretic Complexity–Expressivity Separation", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n}$ be a dataset of size $n >= 4$. The following three sharp complexity separations hold across the 200-algorithm universe:
  1. *(Coordinate-Reading Lower Bound for Prototype Partitioning):* Any deterministic or randomized algorithm achieving a finite multiplicative approximation ratio $alpha < infinity$ to Euclidean $k$-means ($k >= 2$) on arbitrary $cal(X) subset.eq bb(R)^d$ must inspect $Omega(n d)$ scalar coordinate entries in expectation (or pay $Omega(n d)$ preprocessing to construct a coreset).
  2. *(Adversarial Metric-Oracle $Omega(n^2)$ Barrier for Exact Non-Convex Connectivity):* Let $(cal(X), delta)$ be an arbitrary finite metric space accessible only through a symmetric distance oracle returning $delta(x_i, x_j)$ in $O(1)$ time. Any deterministic algorithm that recovers the exact Single Linkage dendrogram (Alg. 56), exact fixed-$epsilon$ DBSCAN connected components with $"minPts" = 2$ (Alg. 31), or the exact 1-medoid / PAM optimum (Alg. 4) must issue at least $binom(n, 2) - floor(n / 2) = Omega(n^2)$ oracle queries in the worst case.
  3. *(Exponential Lattice Barrier for Exhaustive Bottom-Up Subspace Discovery):* Any algorithm required to enumerate all dense axis-parallel subspace clusters (CLIQUE Alg. 134, SUBCLU Alg. 137, ENCLUS Alg. 143) on $cal(X) subset.eq bb(R)^d$ containing a single $s$-dimensional dense cluster with $s <= d$ must examine at least $2^s - 1$ distinct nonempty subspaces in the worst case, imposing an unavoidable $Omega(2^s)$ runtime floor whenever output completeness is enforced.
]

#proof[
  *(1) Coordinate Lower Bound:* Suppose a deterministic algorithm inspects fewer than $(n d) / 4$ coordinates across $cal(X) subset.eq bb(R)^d$ for $k = 2$. Then there exists at least one coordinate pair $(i^*, r^*)$ that is never queried. Consider two datasets $cal(X)_0$ and $cal(X)_M$ identical everywhere except at $(x_(i^*))_(r^*)$, where $cal(X)_0$ places all $n$ points at two distinct clusters ${bold(0), bold(e)_1}$ (so $"OPT"_2 (cal(X)_0) = 0$) whereas $cal(X)_M$ sets $(x_(i^*))_(r^*) = M > 0$ for an arbitrarily large $M$. Because the algorithm never reads $(x_(i^*))_(r^*)$, it outputs the exact same partition $cal(C)$ on both instances. If $cal(C)$ places $x_(i^*)$ in a cluster containing at least one other point, its squared distortion on $cal(X)_M$ is $>= M^2 / 2$, whereas placing $x_(i^*)$ in its own singleton cluster and merging ${bold(0), bold(e)_1}$ yields $"OPT"_2 (cal(X)_M) <= n / 4$. Letting $M -> infinity$ makes the approximation ratio $M^2 / (n / 2) -> infinity$. By Yao's minimax principle applied to the uniform distribution over the uninspected index $(i^*, r^*)$, any randomized algorithm achieving expected ratio $alpha < infinity$ must query $Omega(n d)$ entries.

  *(2) Adversarial Metric-Oracle $Omega(n^2)$ Barrier:* Construct an adversary answering every distance query $delta(x_i, x_j)$ ($i != j$) consistently with the metric axioms by returning $delta(x_i, x_j) = 2$ for every queried pair. Since any assignment of pairwise distances in $[1, 2]$ automatically satisfies the triangle inequality ($1 + 1 >= 2$), the adversary remains completely free to set any *unqueried* pair $(u, v)$ to either $delta_1 (u, v) = 2$ or $delta_2 (u, v) = 1$ without violating any metric axiom.
  - *For Single Linkage (Alg. 56) and DBSCAN with $epsilon = 1.5, "minPts" = 2$ (Alg. 31):* If the algorithm terminates leaving even *one* pair $(u, v)$ unqueried, under metric $delta_1$ all pairs have distance $2 > epsilon$ (producing $n$ isolated singletons and minimum merge height $2$), whereas under $delta_2$ we have $delta_2 (u, v) = 1 <= epsilon$ (merging ${u, v}$ at height $1$). Because the algorithm outputs the exact same result on $delta_1$ and $delta_2$, it errs on at least one valid metric space; hence all $binom(n, 2)$ pairs must be queried.
  - *For exact 1-medoid / PAM (Alg. 4):* Let $U$ be the graph of unqueried pairs upon termination when the algorithm returns candidate medoid $m^*$. Every queried pair has distance $2$. If $U$ contains any edge $(u, v)$ with $u != m^*$ and $v != m^*$ (or any vertex $u != m^*$ incident to more unqueried edges than $m^*$), the adversary sets $delta(u, v) = 1$ and all other unqueried edges to $2$, giving $u$ total distance $2(n - 1) - 1 < 2(n - 1)$ (strictly beating $m^*$ whose distances remain $2$); conversely, if $m^*$ itself has an unqueried edge $(m^*, w)$, setting $delta(m^*, w) = 2$ while setting another vertex's unqueried edge to $1$ again invalidates $m^*$. Thus $U$ cannot contain any edge disjoint from $m^*$ and can contain at most $floor(n / 2)$ matching-bounded edges without creating a strictly superior rival, enforcing at least $binom(n, 2) - floor(n / 2) = Omega(n^2)$ mandatory oracle queries.

  *(3) Exponential Bottom-Up Subspace Barrier:* Let $C subset.eq cal(X)$ be a dense cluster supported on coordinate subset $S subset.eq {1, dots, d}$ with $|S| = s$ and density exceeding threshold $tau$ in a single grid cell (CLIQUE) or $epsilon$-connected ball (SUBCLU). By the Apriori downward-closure monotonicity property (proved in Part VI, Algs. 134 and 137): for *every* nonempty subset $S' subset.eq S$, the coordinate projection $pi_(S') (x)$ satisfies $norm(pi_(S') (x) - pi_(S') (y))_p <= norm(pi_S (x) - pi_S (y))_p$ and maps every cell/ball of $S$ into a containing cell/ball of $S'$ with count $>= |C|$. Therefore $C$ forms a valid dense cluster in *every one* of the $2^s - 1$ nonempty subsets $S' subset.eq S$. Any complete enumerator must output—and hence visit—at least $2^s - 1 = Omega(2^s)$ subspaces. Escaping this barrier strictly requires abandoning exhaustive lattice enumeration in favor of top-down disjoint projection (PROCLUS Alg. 136, ORCLUS Alg. 138), statistical core merging (P3C Alg. 148), or global $ell_1$ self-expression (SSC Alg. 150).
]

Theorem 1 clarifies why *spatial indexing* (Tier 2) and *coresets* (Tier 1) do not contradict the $Omega(n^2)$ barrier: both exploit additional Euclidean or low doubling-dimension geometry unavailable to a general metric oracle—and the moment ambient dimension $d$ grows so that pairwise distances concentrate toward an equidistant simplex @beyer1999nearest, the $Omega(n^2)$ barrier reasserts itself with full force.

#sec("III", "Deterministic vs. Stochastic Trade-Offs and Smoothed Complexity")

A second fundamental dichotomy cutting across Tables IV–XI is *algorithmic determinism vs. stochasticity*. Across the 200 methods, exactly *68 (34.0%)* are strictly deterministic and order-invariant, *31 (15.5%)* are order-dependent deterministic, *65 (32.5%)* are deterministic conditioned on a random seed or subsample, and *36 (18.0%)* evolve via intrinsically stochastic dynamics.

Why do practitioners frequently prefer randomized `Seed` / `Stoch.` algorithms (such as K-Means++, Leiden, Infomap, LDA, and SwAV) over strictly deterministic counterparts (such as Global K-Means, Ward, or Girvan–Newman)? The answer lies in a deep mathematical tension between *reproducibility*, *local-minimum escape*, and *worst-case vs. smoothed runtime* @arthur2007kmeanspp, @arthur2011smoothed, @vattani2011kmeans.

#theorem(title: "Determinism–Approximation Trade-Off, Multi-Start Concentration & Smoothed Complexity", origin: "Original Unified Proof")[
  Let $cal(X) subset.eq bb(R)^d$ with $|cal(X)| = n$ and target cluster count $k >= 2$, and denote the global $k$-means optimum by $"OPT"_k = min_(C : |C| = k) Phi(cal(X); C)$. Then:
  1. *(Deterministic Local-Optimum Trap vs. Quadratic Global Premium):* Deterministic single-run Lloyd refinement (Alg. 1) can terminate at local minima with arbitrarily bad ratio $Phi / "OPT"_k -> infinity$ and can require $T = 2^(Omega(sqrt(n)))$ iterations in the worst case @vattani2011kmeans. Eliminating seeding randomness deterministically via Global K-Means (Alg. 17 @likas2003global) escalates per-dataset cost from $O(n k d T)$ to $O(n^2 k^2 d T)$—a multiplicative factor of $Theta(n k)$.
  2. *(Exponential High-Probability Concentration of Multi-Start $D^2$-Seeding):* Executing $R$ independent runs of K-Means++ seeding (Alg. 2 @arthur2007kmeanspp) and returning the minimum-cost solution $hat(C)_R$ achieves, for any failure tolerance $delta in (0, 1)$,
     $ bb(P)(Phi(cal(X); hat(C)_R) <= 16 (ln k + 2) "OPT"_k) >= 1 - delta $
     whenever $R >= ceil(log_2 (1 / delta))$, in total wall time $O(R n k d T) = O(n k d T log(1 / delta))$—strictly beating Global K-Means by a factor of $Theta(n k / log(1 / delta))$.
  3. *(Smoothed Polynomial Reconciliation):* Under Gaussian coordinate perturbation $tilde(x)_i = x_i + g_i$ with $g_i tilde.op cal(N)(bold(0), sigma^2 I_d)$ and bounded diameter $norm(x_i)_2 <= 1$, the expected number of Lloyd iterations satisfies the Arthur–Manthey–Röglin smoothed polynomial bound @arthur2011smoothed
     $ bb(E)[T_"Lloyd" (tilde(cal(X)))] = O("poly"(n, d, 1 / sigma)), $
     reconciling Vattani's exponential worst-case deterministic lower bound with empirically bounded $T approx 10$–$100$ on real noisy data.
]

#proof[
  *(1) Deterministic Trap and Premium:* Place $k = 2$ clusters in $bb(R)$ at ${0, 0}$ and ${M, M + 1}$ with $M >> 1$ (so $"OPT"_2 = 1 / 2$ achieved by centers ${0, M + 1 / 2}$). If a deterministic initialization places one center at $M$ and the other at $M + 1$, the Voronoi partition splits ${M + 1}$ from ${0, 0, M}$ with stationary Lloyd centers ${M / 3, M + 1}$ and cost $Phi = (2 / 3) M^2 -> infinity$ as $M -> infinity$, yielding an unbounded approximation ratio. Vattani @vattani2011kmeans constructs an explicit planar $d = 2$ gadget of $n$ points where exact Lloyd iterations visit $2^(Omega(sqrt(n)))$ distinct Voronoi configurations. To remove seeding variance deterministically, Global K-Means (Alg. 17) adds centers sequentially $m = 1, dots, k$: at stage $m$, it tests every point $x_i in cal(X)$ ($n$ candidates) as the $m$-th initial centroid and runs full Lloyd refinement ($O(n m d T)$ per candidate), summing to
  $ sum_(m=1)^k n dot.op O(n m d T) = O(n^2 k^2 d T). $

  *(2) Multi-Start $D^2$-Seeding Concentration:* By the Arthur–Vassilvitskii theorem @arthur2007kmeanspp (proved in Part I, Alg. 2), a single independent draw $C^((r))$ of $D^2$-weighted seeding satisfies
  $ bb(E)[Phi(cal(X); C^((r)))] <= 8 (ln k + 2) "OPT"_k. $
  Because $Phi(cal(X); C^((r))) >= "OPT"_k > 0$ is a nonnegative random variable, Markov's inequality yields
  $ bb(P)(Phi(cal(X); C^((r))) > 16 (ln k + 2) "OPT"_k) \
    <= (bb(E)[Phi(cal(X); C^((r)))]) / (16 (ln k + 2) "OPT"_k) <= 1 / 2. $
  Let $hat(C)_R := arg min_(1 <= r <= R) Phi(cal(X); C^((r)))$, optionally refined by Lloyd's algorithm (which monotonically decreases $Phi$, so $Phi_"refined" <= Phi(cal(X); hat(C)_R)$). By independence of the $R$ runs,
  $ bb(P)(Phi(cal(X); hat(C)_R) > 16 (ln k + 2) "OPT"_k) \
    = product_(r=1)^R bb(P)(Phi(cal(X); C^((r))) > 16 (ln k + 2) "OPT"_k) <= 2^(-R). $
  Choosing $R >= ceil(log_2 (1 / delta))$ enforces $2^(-R) <= delta$, proving the stated high-probability guarantee in $O(R n k d T) = O(n k d T log(1 / delta))$ time (e.g., $R = 20$ drives failure probability below $10^(-6)$).

  *(3) Smoothed Polynomial Reconciliation:* In Arthur–Manthey–Röglin @arthur2011smoothed, every Lloyd iteration reduces the discrete potential $Phi$ by at least $Delta_t = sum_j |C_j^((t))| norm(mu_j^((t)) - mu_j^((t-1)))_2^2$. The total potential starts bounded by $O(n (1 + d sigma^2))$ with high probability. The adversarial exponential chain in @vattani2011kmeans requires exponentially tiny potential drops $Delta_t = 2^(-Omega(sqrt(n)))$ realized by delicate co-linear alignments. Under i.i.d. Gaussian noise $cal(N)(bold(0), sigma^2 I_d)$, the probability that any sequence of $O(k)$ consecutive swaps alters the centroid geometry while dropping $Phi$ by less than $epsilon_0 = 1 / "poly"(n, d, 1 / sigma)$ is bounded via anti-concentration of linear combinations of Gaussians, bounding the expected number of iterations by $O(n^(34) k^(34) d^8 log^4 (n) / sigma^6)$—strictly polynomial in $(n, d, 1 / sigma)$.
]

#subsec("A", "The Hidden Fragility of Order-Dependent Deterministic Algorithms")

Thirty-one algorithms across Tables IV–XI fall into the subtle `Order` class: they invoke *zero* pseudo-random number generator calls—giving the illusion of full determinism—yet their output depends sensitively on the arrival permutation $pi in S_n$ of the input records. Prime examples span all families: MacQueen's online K-Means (Alg. 21), COP-KMeans (Alg. 20), DBSCAN border-point ties (Alg. 31), BIRCH CF-tree insertion (Alg. 64), BUBBLE (Alg. 71), PERCH/GRINCH (Algs. 76–77), COBWEB/CLASSIT (Algs. 85–86), Louvain node-visit order (Alg. 107), ART vigilance ordering (Alg. 178), and streaming engines CluStream, DenStream, ClusTree, and D-Stream (Algs. 193–194, 196–197). The next proposition proves an explicit permutation-sensitivity bound—and shows how coreset merge-and-reduce eliminates it.

#proposition(title: "Order-Sensitivity Gap & Coreset Permutation Stabilization", origin: "Original Proof")[
  1. *(Unbounded Distortion & Cardinality Drift of Fixed-Threshold Incremental Trees):* Consider single-pass diameter-bounded BIRCH (Alg. 64 @zhang1996birch) with diameter threshold $D_"max" > 0$ (equivalently ART Alg. 178 @carpenter1987art with $L_1$ box width $D_"max"$). For every integer $m >= 1$ and every $delta in (0, 1 / 10)$, there exists a 1D dataset $cal(X) subset.eq bb(R)$ of $n = 4 m$ points and two arrival permutations $pi_1, pi_2 in S_n$ such that:
     $ |cal(L)(pi_2)| = 3 / 2 |cal(L)(pi_1)| = 3 m, \
       ("SSE"(pi_2)) / ("SSE"(pi_1)) = (2 (1 - delta)^2) / delta^2 -> infinity quad (delta -> 0^+). $
  2. *(Coreset Permutation Stabilization):* Conversely, let $cal(A)$ be a merge-and-reduce $(epsilon, k)$-coreset streaming algorithm (BICO Alg. 78 @fichtenberger2013bico, StreamKM++ Alg. 195 @ackermann2012streamkmpp) producing weighted coreset $(S_pi, w_pi)$ satisfying $(1 - epsilon) Phi(cal(X); C) <= Phi(S_pi, w_pi; C) <= (1 + epsilon) Phi(cal(X); C)$ for all $k$-center sets $C$ and every arrival permutation $pi in S_n$. Then for *any two permutations* $pi_1, pi_2 in S_n$, exact offline $k$-means solutions $C_1^*, C_2^*$ computed on $S_(pi_1)$ and $S_(pi_2)$ satisfy
     $ (Phi(cal(X); C_1^*)) / (Phi(cal(X); C_2^*)) <= (1 + epsilon) / (1 - epsilon) \
       <= ((1 + epsilon) / (1 - epsilon))^2 = 1 + 4 epsilon + O(epsilon^2). $
]

#proof[
  *(1) Constructive Permutation Counterexample:* Set threshold $D_"max" = 2 - (3 / 2) delta$ with $delta in (0, 1 / 10)$, and construct $m$ isolated four-point blocks $B_j = {10 j, 10 j + delta, 10 j + 2 - delta, 10 j + 2}$ for $j = 0, dots, m - 1$. Because distinct blocks are separated by $>= 8 > D_"max"$, points from different blocks can never merge into the same leaf CF-entry. Inside a single block $B_0 = {0, delta, 2 - delta, 2}$, observe the exact pairwise diameters:
  - $"diam"({0, delta}) = "diam"({2 - delta, 2}) = delta <= D_"max"$ (valid merges, each with SSE $= delta^2 / 2$);
  - $"diam"({delta, 2 - delta}) = 2 - 2 delta < 2 - (3 / 2) delta = D_"max"$ (valid merge with SSE $= 2 (1 - delta)^2$);
  - $"diam"({0, delta, 2 - delta}) = "diam"({delta, 2 - delta, 2}) = 2 - delta > D_"max"$ and $"diam"({0, 2}) = 2 > D_"max"$ (all invalid).

  Now compare two presentation orders applied block-by-block:
  - *Order $pi_1 = (0, delta, 2 - delta, 2)$:* Point $0$ opens leaf $L_1 = {0}$; point $delta$ merges into $L_1$ ($delta <= D_"max"$); point $2 - delta$ cannot join $L_1$ ($2 - delta > D_"max"$) and opens $L_2 = {2 - delta}$; point $2$ joins its nearest leaf $L_2$ ($delta <= D_"max"$). Thus each block yields $2$ leaves with block SSE $= delta^2 / 2 + delta^2 / 2 = delta^2$, giving $|cal(L)(pi_1)| = 2 m$ and $"SSE"(pi_1) = m delta^2$.
  - *Order $pi_2 = (delta, 2 - delta, 0, 2)$:* Point $delta$ opens $L_1 = {delta}$; point $2 - delta$ merges into $L_1$ since $2 - 2 delta <= D_"max"$, forming $L_1 = {delta, 2 - delta}$; point $0$ cannot merge into $L_1$ ($2 - delta > D_"max"$) and opens singleton $L_2 = {0}$; point $2$ cannot merge into $L_1$ ($2 - delta > D_"max"$) nor $L_2$ ($2 > D_"max"$) and opens singleton $L_3 = {2}$. Thus each block yields $3$ leaves with block SSE $= 2 (1 - delta)^2 + 0 + 0$, giving $|cal(L)(pi_2)| = 3 m = (3 / 2) |cal(L)(pi_1)|$ and $"SSE"(pi_2) / "SSE"(pi_1) = 2 (1 - delta)^2 / delta^2 -> infinity$ as $delta -> 0^+$.

  *(2) Coreset Permutation Stability:* Let $C_"opt" = arg min_C Phi(cal(X); C)$. For any permutation $pi_r$ ($r in {1, 2}$), the $(epsilon, k)$-coreset inequalities give
  $ (1 - epsilon) Phi(cal(X); C_r^*) &<= Phi(S_(pi_r), w_(pi_r); C_r^*) \
    &<= Phi(S_(pi_r), w_(pi_r); C_"opt") <= (1 + epsilon) Phi(cal(X); C_"opt"). $
  Dividing by $1 - epsilon$ yields $"OPT"_k <= Phi(cal(X); C_r^*) <= ((1 + epsilon) / (1 - epsilon)) "OPT"_k$ for both $r = 1, 2$. Taking the ratio between $r = 1$ and $r = 2$ bounds $Phi(cal(X); C_1^*) / Phi(cal(X); C_2^*) <= (1 + epsilon) / (1 - epsilon)$ (and even with $(1 + epsilon)$-approximate offline solvers on each coreset, the ratio remains strictly bounded by $((1 + epsilon) / (1 - epsilon))^2 = 1 + 4 epsilon + O(epsilon^2)$).
]

Proposition 1 delivers an immediate engineering rule: *never deploy an order-dependent `Order` algorithm on pre-sorted database tables without either random shuffling (converting `Order` into low-variance `Seed`) or coreset stabilization (BICO / #box[StreamKM++])*.

#sec("IV", "Memory Footprints, Working-Set Hierarchies & High-Dimensional Scaling")

As dataset sizes have grown from $n = 10^3$ to $n = 10^8$, *working-set RAM*—not asymptotic FLOP count—has become the decisive failure boundary in production clustering systems. Stratifying all 200 methods in Tables IV–XI by peak working memory yields five cleanly separated footprint classes, cross-tabulated against determinism class in Table III.

#wide_figure(caption: [*Table III:* Joint Memory-Footprint vs. Determinism-Class Contingency Matrix Across Algorithms 1–200, with Practical Single-Node Scale Ceilings (64 GB RAM, IEEE-754 float64). Every cell count reconciles strictly to Tables IV–XI.])[
  #set par(justify: false, first-line-indent: 0pt)
  #set text(size: 7.1pt)
  #show math.equation: set text(size: 7.0pt)
  #table(
    columns: (1.55fr, 1.15fr, 0.42fr, 0.42fr, 0.44fr, 0.46fr, 0.56fr, 1.25fr),
    align: (left, left, center, center, center, center, center, left),
    inset: (x: 3.4pt, y: 2.9pt),
    fill: (_, y) => if y == 0 { rgb("#f1f5f9") } else if y == 6 { rgb("#eff6ff") } else if calc.odd(y) { rgb("#fbfdff") } else { white },
    table.header(
      [Memory-Footprint Class], [Working-Set Scaling], [Det.], [Order], [Seed], [Stoch.], [Total (%)], [64 GB Single-Node Scale],
    ),
    [*Sublinear / Constant-in-$n$ Streaming*], [$O(k d)$, $O(q d)$, $O(S d)$, $O(K)$, $O(s^2)$], [2], [13], [5], [2], [*22 (11.0%)*], [Unbounded ($n -> infinity$, stream)],
    [*Linear Coordinate / Sparse CSR Graph*], [$O(n d)$, $O(n + m)$, $O(N)$], [26], [15], [24], [13], [*78 (39.0%)*], [$n approx 10^7$–$10^9$ (RAM / CSR)],
    [*Bilinear Responsibility / Two-Mode*], [$O(n (d + k))$, $O(n p)$, $O(k (n + V))$], [11], [2], [24], [7], [*44 (22.0%)*], [$n approx 10^6$–$10^7$ (moderate $k$)],
    [*Quadratic Pairwise / Gram / Co-Assoc.*], [$O(n^2)$, $O(V n^2)$], [28], [1], [4], [1], [*34 (17.0%)*], [$n approx 10^4$–$1.4 times 10^5$ (RAM wall)],
    [*Parametric Covariance / Deep Network*], [$O(n k + k d^2)$, $O(|theta| + B r)$], [1], [0], [8], [13], [*22 (11.0%)*], [$n approx 10^6$–$10^8$ (GPU / small $d$)],
    [*Total Across All 200 Methods*], [—], [*68*], [*31*], [*65*], [*36*], [*200 (100%)*], [*Span: $10^4$ to unbounded stream*],
  )
]

Reading Table III vertically exposes a profound systems asymmetry:
- *Sublinear Streaming Engines (22 algs., 11.0%):* Maintain working state strictly independent of $n$ (or $O(log n)$ coreset buckets), enabling infinite-horizon deployment on $n > 10^9$ streams—yet *only two* strictly deterministic (`Det.`) methods achieve $o(n)$ RAM (fixed spatial grids STING Alg. 133 and WaveCluster Alg. 135).
- *Linear Coordinate & Sparse-CSR Engines (78 algs., 39.0%):* Represent the largest single footprint tier, holding the raw data $O(n d)$ or sparse graph CSR arrays $O(n + m)$ without ever materializing an $n times n$ affinity matrix.
- *Bilinear & Two-Mode Matrices (44 algs., 22.0%):* Store soft responsibility tables $U in [0, 1]^(n times k)$ (FCM, PCM, LCA, Poisson/Bernoulli mixtures, topic models, biclustering matrices) or $k$-NN lists ($O(n k)$). Memory remains linear in $n$ for fixed $k$, but blows up under extreme clustering ($k approx sqrt(n)$), motivating sparse top-$r$ responsibility pruning or online micro-batches.
- *The Quadratic $O(n^2)$ RAM Wall (34 algs., 17.0%):* Thirty-four algorithms—including PAM/FasterPAM, Kernel K-Means, KFCM, dense UPGMA/Ward/DIANA/Minimax/NJ, dense Spectral NJW/Ncut/Ratio Cut, Affinity Propagation, DPC, Walktrap, SSC, Monti Consensus, and CSPA—require $Theta(n^2)$ pairwise storage. Notably, *28 / 34 (82.4%)* of these quadratic-memory methods are strictly deterministic (`Det.`), confirming that exact global pairwise determinism carries a severe memory tax.
- *Parametric Covariance & Deep Neural State (22 algs., 11.0%):* Splits cleanly into two opposite scaling behaviors: full-covariance mixtures (GMM, DPGMM, TMM, MCLUST, Mixmod, Skew-Normal, IGMM, GK, GG, ORCLUS) pay $O(k d^2)$ memory that explodes in high dimension $d$, whereas mini-batch deep representation clusterers (DEC, IDEC, DCN, VaDE, DEPICT, DAC, IIC, SwAV, CC) pay $O(|theta| + B d)$ VRAM—decoupled from dataset size $n$ (unlike offline feature-caching DeepCluster Alg. 191, SCAN Alg. 189, and JULE Alg. 186, which store $O(n r)$ embeddings or $O(n k)$ neighbor indices).

The next theorem formalizes the exact byte-level RAM wall confronting the 34 quadratic-memory algorithms and proves three principled mathematical escapes.

#theorem(title: "The Quadratic RAM Wall & Coreset / Sparse-Graph Escape Bounds", origin: "Original Unified Proof")[
  Let $cal(X) = {x_1, dots, x_n}$ be processed on a single machine with usable RAM budget $B_"RAM"$ bytes (IEEE-754 double precision, $8$ bytes per scalar). Then:
  1. *(Exact $O(n^2)$ Memory Ceiling):* Any algorithm materializing a symmetric zero-diagonal pairwise dissimilarity, kernel Gram, or co-association matrix $M in bb(R)^(n times n)$ in compact upper-triangular float64 storage requires at least $M_"sym" (n) = 4 n (n - 1)$ bytes (and $8 n^2$ bytes for full dense BLAS/LAPACK eigensolvers). Consequently, the hard feasibility ceiling satisfies
     $ n_"max" (B_"RAM") <= floor(1 / 2 + sqrt(1 / 4 + B_"RAM" / 4)) \
       approx 0.5 sqrt(B_"RAM") quad ("compact"), $
     capping a 64 GB workstation ($B_"RAM" = 6.4 times 10^(10)$ bytes) at $n <= 1.26 times 10^5$ (compact) or $n <= 8.94 times 10^4$ (full dense eigensolve).
  2. *(Coreset Memory Collapse for Prototype & Variance Objectives):* For Euclidean $k$-means (Algs. 1–2, 59, 64, 78, 195), Har-Peled–Mazumdar @harpeled2004coresets and BICO @fichtenberger2013bico construct an $(epsilon, k)$-coreset of cardinality $S = O(k epsilon^(-d) log n)$ (or dimension-independent $S = O(k d epsilon^(-2) log k)$ via sensitivity sampling), compressing RAM from $Theta(n d)$ (or $Theta(n^2)$ for Ward) down to $O(S d) = o(n)$ while guaranteeing $(1 + epsilon)$-optimal distortion.
  3. *(Sparse $k$-NN CSR & Nyström Escape for Spectral / Kernel Methods):* Replacing a dense $n times n$ Gaussian affinity matrix $W$ by a symmetric $m_"nn"$-nearest-neighbor CSR graph reduces memory from $8 n^2$ bytes to $approx 16 n m_"nn" + 8 n$ bytes—a compression factor of $n / (2 m_"nn")$ (exceeding $2,500 times$ at $n = 10^5, m_"nn" = 20$)—while rank-$r$ Nyström approximation $hat(K) = C W_r^(+) C^T$ cuts Kernel K-Means / Spectral memory from $O(n^2)$ to $O(n r + r^2)$ with operator-norm error bounded by $norm(K - hat(K))_2 <= lambda_(r + 1) (K) + O(norm(K)_F / sqrt(r))$.
]

#proof[
  *(1) Exact Byte Ceiling:* A symmetric zero-diagonal matrix $M in bb(R)^(n times n)$ contains exactly $binom(n, 2) = n (n - 1) / 2$ independent floating-point entries. At $8$ bytes per IEEE-754 binary64 scalar, storing the packed triangle demands $8 dot.op n (n - 1) / 2 = 4 n (n - 1)$ bytes. Solving $4 n^2 - 4 n <= B_"RAM"$ via the quadratic formula yields $n <= (1 + sqrt(1 + B_"RAM")) / 2$. Substituting $B_"RAM" = 6.4 times 10^(10)$ bytes gives $n <= 126,491$ for packed symmetric matrices (PAM, SLINK/UPGMA matrix mode, DPC, CSPA), while symmetric dense eigensolvers (NJW, Ncut, Ratio Cut, SSC) require a full $n times n$ matrix ($8 n^2$ bytes plus LAPACK workspace), capping $n <= sqrt(B_"RAM" / 8) approx 89,442$. At $n = 10^6$, even packed symmetric storage requires $4.0$ TB of RAM—completely exceeding commodity nodes.

  *(2) Coreset Compression:* Let $(S, w)$ be an $(epsilon, k)$-coreset of $cal(X)$ satisfying $|Phi(S, w; C) - Phi(cal(X); C)| <= epsilon Phi(cal(X); C)$ uniformly over all $k$-center configurations $C subset.eq bb(R)^d$. By the additive Bregman sufficient-statistic identity (Comparison Chapter I, Theorem 1), any streaming merge-and-reduce binary tree of depth $L = ceil(log_2 (n / m))$ combining buckets of size $m$ with per-level accuracy $epsilon_L = epsilon / (2 L)$ requires storing at most $L + 1$ active buckets simultaneously—yielding total working RAM $O(m d log(n / m))$ (StreamKM++ Alg. 195 @ackermann2012streamkmpp) or fixed-budget $O(S d)$ via adaptive threshold doubling (BICO Alg. 78 @fichtenberger2013bico)—both strictly $o(n)$ as $n -> infinity$, while maintaining total relative distortion error $(1 + epsilon / (2 L))^L <= e^(epsilon / 2) <= 1 + epsilon$ for $epsilon in (0, 1)$.

  *(3) Sparse CSR and Nyström Compression:* Storing an undirected $m_"nn"$-NN graph with at most $n m_"nn"$ undirected edges in Compressed Sparse Row (CSR) format requires $2 n m_"nn"$ column indices ($4$ bytes each in `int32`), $2 n m_"nn"$ edge weights ($8$ bytes in `float64`, or $4$ bytes in `float32`), and $n + 1$ row pointers ($8$ bytes in `int64`), totaling $<= 24 n m_"nn" + 8 n$ bytes (or $16 n m_"nn" + 8 n$ bytes in `float32` weights). Dividing $8 n^2$ by $16 n m_"nn"$ yields compression ratio $n / (2 m_"nn")$ (at $n = 10^5, m_"nn" = 20$, $80$ GB shrinks to $32$ MB). For Nyström approximation, sampling $r << n$ landmark columns yields $C in bb(R)^(n times r)$ and $W_r in bb(R)^(r times r)$, storing $O(n r + r^2)$ scalars; the spectral-norm residual decomposes into the optimal rank-$r$ tail $lambda_(r + 1) (K)$ plus the standard leverage-score sampling perturbation $O(norm(K)_F / sqrt(r))$, reducing eigendecomposition from $O(n^3)$ to $O(n r^2 + r^3)$.
]

#corollary(title: "Ten-Dimension Pareto Frontier & No-Free-Lunch Conservation", origin: "Original Proof")[
  Across all 200 algorithms documented in Tables IV–XI, *no single algorithm simultaneously achieves*: (i) sublinear $o(n)$ working memory, (ii) strict permutation-invariant determinism (`Det.`), (iii) arbitrary non-convex manifold recovery, (iv) explicit breakdown robustness strictly exceeding $1 / n$, (v) automatic model-order $k$ selection, and (vi) strong high-dimensional survival as $d -> infinity$. Every method in Parts I–VIII lies on a constrained trade-off surface governed by three conservation laws:
  1. *Memory vs. Exact Non-Convex Determinism:* Achieving (ii) and (iii) on arbitrary metrics forces $Omega(n^2)$ queries (Theorem 1(2)) and $Omega(n)$–$Omega(n^2)$ state (violating (i)).
  2. *Streaming Additivity vs. Outlier Breakdown:* Achieving (i) via additive sufficient statistics (BIRCH, BICO, CluStream) locks the breakdown point to $1 / n$ (violating (iv)), whereas bounded-influence trimming/medoids (PAM, TMM, Trimmed K-Means) break additive sufficiency.
  3. *High-$d$ Concentration vs. Nonparametric Level Sets:* Achieving (iii) and (v) via Euclidean density estimation (Cat. II) suffers exponential sample-complexity decay $O(n^(-4 / (d + 4)))$ (violating (vi)), forcing high-$d$ survival into directional normalization (Cat. I/IV), sparse graphs (Cat. V), subspace projection (Cat. VI), or learned neural encoders (Cat. VIII).
]

#proof[
  Inspecting the 200 rows of Tables IV–XI directly—or combining Theorem 1, Proposition 1, Theorem 3, Kleinberg's Impossibility Theorem @kleinberg2002impossibility, and Beyer's distance-concentration theorem @beyer1999nearest—verifies all three exclusions simultaneously. Specifically: (1) Only two strictly deterministic methods sit in the sublinear-RAM class (Table III: STING Alg. 133 and WaveCluster Alg. 135), and both rely on fixed $c^d$ grids that fail exponentially in high $d$ (violating (vi)). (2) Every algorithm with Strong high-$d$ resilience and arbitrary non-convex recovery either requires a sparse graph ($Omega(n + m)$ RAM, Cat. V), subspace search ($Omega(n d)$ to $O(n^2)$ RAM, Cat. VI), or stochastic deep neural training (`Stoch.`, violating (ii) and (v), Cat. VIII). (3) Every strictly deterministic `Det.` algorithm with Auto-$k$ and Strong high-$d$ resilience (SNN Alg. 36, SNN-DPC Alg. 52, SCC Alg. 79, FINCH Alg. 80, COSA Alg. 145, FIRES Alg. 146, INSCY Alg. 147) requires at least $Omega(n d)$ to $Omega(n^2)$ memory (violating (i)). Hence the six-property intersection is mathematically empty across all 200 algorithms.
]

#subsec("A", "Cache-Hierarchy Locality, High-Dimensional Phase Transitions & Lookup Protocol")

Beyond asymptotic RAM bounds, modern CPU/GPU execution speed across Tables IV–XI is governed by *memory-access regularity* and *dimension-conditional index collapse*:

1. *Sequential BLAS-3 vs. Pointer-Chasing Cache Misses:* Bilinear Tier 3 algorithms operating on contiguous dense matrices—K-Means (Alg. 1), Spherical K-Means (Alg. 12), GMM / Poisson / Bernoulli EM (Algs. 81, 95–96), NMF (Alg. 92), FCM (Alg. 151), and mini-batch deep encoders (Algs. 181–192)—cast inner loops as Level-3 GEMM matrix multiplications $X M^T$, sustaining $> 80%$ of peak SIMD/tensor-core FLOPs with minimal L3 cache misses. By contrast, nominally faster Tier 1–2 sparse-graph and spatial-tree algorithms—Louvain (Alg. 107), Leiden (Alg. 108), LPA (Alg. 110), $k$-d tree DBSCAN (Alg. 31), and ClusTree (Alg. 196)—execute irregular indirect memory fetches (`adj[ptr[v]]`) that stall on DRAM latency ($approx 60$–$100$ ns per miss). Consequently, for small-to-medium $n <= 5 times 10^4$, cache-resident Tier 3 vector kernels routinely outrun pointer-chasing Tier 2 trees by $3 times$–$10 times$ in wall-clock time—an empirical crossover directly verified in Comparison Chapter XIII.
2. *The $d_"crit" approx 15$–$20$ Spatial-Index Collapse Boundary:* Every algorithm marked Tier 2 ($O(n log n)$) with `Poor` high-$d$ rating in Tables V and IX relies on axis-aligned bounding-box pruning ($k$-d tree, $R^*$-tree, quadtree, or regular cell stencil). As proved via Beyer's concentration theorem @beyer1999nearest, once intrinsic dimension exceeds $d_"crit" approx 15$–$20$, an $epsilon$-query sphere intersects almost every leaf bounding box in the tree, forcing $Theta(n)$ leaf visits per point and collapsing effective runtime from $O(n log n)$ straight to $O(n^2 d)$ with a worse constant factor than a SIMD brute-force distance scan @schubert2017dbscanrevisited.
3. *Four-Step Diagnostic Protocol for Navigating Tables IV–XI:* When consulting Section V to shortlist candidate algorithms for a new dataset $(cal(X), d, n)$, filter sequentially along four orthogonal columns: *(i) Space & Time Scale*—eliminate any method whose Space column hits Theorem 3's $O(n^2)$ RAM wall when $n > 10^5$, or whose Time column is Tier 5–6; *(ii) High-$d$ Survival*—when $d > 25$ (without low-dimensional manifold structure), restrict immediately to `Good` or `Strong` rows (directional, sparse-graph, subspace, or deep-encoder methods); *(iii) Geometry & Noise*—match convex vs. non-convex / graph / temporal inductive bias and require explicit noise excision or bounded influence whenever contamination exceeds $1%$; and *(iv) Determinism & Auto-$k$*—choose `Det.` or multi-start `Seed` (Theorem 2) for regulated reproducible pipelines, and enforce `Yes` / `Dendro.` in Auto-$k$ whenever cluster count is unknown a priori.

// =============================================================================
// FULL-WIDTH SINGLE-COLUMN MASTER REFERENCE COMPENDIUM (TABLES IV–XI)
// =============================================================================
#set page(columns: 1)

#sec("V", "Exhaustive 200-Algorithm Master Comparative Reference Matrices (Tables IV–XI)")

Tables IV through XI present the complete, zero-gap 10-dimension reference matrices for every clustering algorithm analyzed in Parts I–VIII (IDs 1–200), organized strictly by the monograph's eight taxonomic categories. Within every table, column headers repeat automatically across page breaks, asymptotic complexity matches the formal proofs in each algorithm's dedicated chapter, and the four determinism codes (`Det.`, `Order`, `Seed`, `Stoch.`) and four high-dimensional tiers (`Poor`, `Fair`, `Good`, `Strong`) follow the standardized definitions of Section I-A.

#subsec("A", "Table IV — Cat. I: Centroid, Medoid & Partitioning (Algorithms 1–30)")

Category I spans thirty prototype-based partitioners anchored by Lloyd's Bregman/Voronoi engine (Alg. 1 @lloyd1982least). Twenty-one methods reside in the $O(n k d T)$ bilinear tier, while Mini-Batch K-Means (Alg. 3 @sculley2010webscale), CLARA / CLARA-MED (Algs. 5, 25 @kaufman1990finding), and parallel K-Means|| (Alg. 23 @bahmani2012scalable) break into linear or sublinear streaming complexity. Only three methods select $k$ automatically (X-Means @pelleg2000xmeans, G-Means @hamerly2003gmeans, ISODATA @ball1965isodata), and raising the $1 / n$ breakdown point requires medoids (Algs. 4–6, 24–25 @schubert2019fasterpam), coordinate $L_1$ medians (Alg. 15 @bradley1997clustering), or impartial $alpha$-trimming (Alg. 26 @cuesta1997trimmed). In high dimensions, Spherical K-Means (Alg. 12 @dhillon2001concept), Power K-Means (Alg. 27 @xu2019power), and feature-weighted $L_p$ variants (Algs. 28–29 @huang2005automated, @deamorim2012minkowski) form the empirical Pareto frontier.

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table IV: Complete 10-Dimension Theoretical Reference Matrix — Cat. I: Centroid, Medoid & Partitioning (Algs. 1–30; $N = 30$, Auto-$k$: 3/30, High-$d$ Good/Strong: 6/30, Strictly Det.: 4/30).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*1*], [K-Means (Lloyd, 1982)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Convex $L_2$ Voronoi], [Fragile ($1 / n$)], [No], [Poor ($L_2$)], [General tabular, VQ],
    [*2*], [K-Means++ (2007)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Convex $L_2$ Voronoi ($D^2$)], [Fragile ($1 / n$)], [No], [Poor ($L_2$)], [Default $k$-means seeding],
    [*3*], [Mini-Batch K-Means (2010)], [$O(b k d T)$], [$O((b + k) d)$], [Stoch.], [Convex $L_2$ Voronoi], [Fragile (SGD avg.)], [No], [Poor ($L_2$)], [Web-scale streaming VQ],
    [*4*], [K-Medoids / PAM (1990)], [$O(k (n - k)^2 T)$], [$O(n^2)$], [Det.], [Arbitrary metric balls], [Robust (medoids)], [No], [Fair (metric)], [Arbitrary metric / bioinfo.],
    [*5*], [CLARA (1990)], [$O(R (k s^2 T + n k))$], [$O(s^2)$ stream], [Seed], [Subsampled metric balls], [Robust (medoids)], [No], [Fair (metric)], [Large-scale medoid data],
    [*6*], [CLARANS (2002)], [$O(r m_"max" n T)$], [$O(n d)$], [Stoch.], [Johnson-graph medoids], [Robust (medoids)], [No], [Fair (metric)], [Spatial DBs, GIS routing],
    [*7*], [K-Modes (1998)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Hamming discrete modes], [Moderate (modes)], [No], [Good (Ham.)], [Pure categorical / survey],
    [*8*], [K-Prototypes (1998)], [$O(n k (d_r + d_c) T)$], [$O((n + k) d)$], [Seed], [Mixed $L_2$ + Hamming], [Sensitive ($L_2$ part)], [No], [Fair (mixed)], [Enterprise CRM / EHR],
    [*9*], [Bisecting K-Means (2000)], [$O(n d R T log k)$], [$O((n + k) d)$], [Seed], [Hierarchical $L_2$ splits], [Fragile ($1 / n$)], [No], [Fair (text)], [Document hierarchies],
    [*10*], [X-Means (2000)], [$O(n k_"max" d T)$], [$O((n + k_"max") d)$], [Seed], [Spherical Gaussian BIC], [Fragile ($1 / n$)], [Yes (BIC)], [Poor (iso.)], [Exploratory tabular DBs],
    [*11*], [G-Means (2003)], [$O(n k (d T + log n))$], [$O((n + k) d)$], [Seed], [1D PCA Gaussian splits], [Sensitive (tails)], [Yes (AD)], [Fair (PCA)], [Gaussian exploratory data],
    [*12*], [Spherical K-Means (2001)], [$O("nnz"(X) k T)$], [$O("nnz"(X) + k d)$], [Seed], [Hypersphere $bb(S)^(d-1)$ cos], [Moderate (normed)], [No], [Strong (cos)], [TF-IDF, dense embeddings],
    [*13*], [Elkan's K-Means (2003)], [$O(n k d T)$ (pruned)], [$O(n k + n d)$], [Seed], [Exact $L_2$ Voronoi], [Fragile ($1 / n$)], [No], [Fair (bounds)], [Exact high-$k$ acceleration],
    [*14*], [Hamerly's K-Means (2010)], [$O(n k d T)$ (pruned)], [$O(n d + k d)$], [Seed], [Exact $L_2$ Voronoi], [Fragile ($1 / n$)], [No], [Poor (low-$d$)], [Low-$d$ memory-tight VQ],
    [*15*], [K-Medians (1997)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Axis-aligned $L_1$ diamonds], [Robust (50% axis)], [No], [Fair ($L_1$)], [Heavy-tailed / finance],
    [*16*], [Kernel K-Means (1998/2004)], [$O(n^2 (d + T))$], [$O(n^2)$], [Seed], [Non-convex RKHS hulls], [Kernel-dependent], [No], [Good (RKHS)], [Manifolds, graph kernels],
    [*17*], [Global K-Means (2003)], [$O(n^2 k^2 d T)$], [$O((n + k) d)$], [Det.], [Deterministic $L_2$ Voronoi], [Fragile ($1 / n$)], [No], [Poor ($L_2$)], [Small reproducible bench.],
    [*18*], [ISODATA (1965)], [$O((n k + k^2) d T)$], [$O((n + k) d + k^2)$], [Seed], [Adaptive $L_2$ split/merge], [Prunes tiny clus.], [Yes (heur.)], [Poor ($L_2$)], [Remote sensing, multispec.],
    [*19*], [Hartigan-Wong (1979)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Strict transfer $L_2$ opt.], [Fragile ($1 / n$)], [No], [Poor ($L_2$)], [R default / fine $L_2$ opt.],
    [*20*], [COP-KMeans (2001)], [$O(T k (n d + |E|))$], [$O(n d + |E|)$], [Order], [Constrained $L_2$ Voronoi], [Fragile + ML/CL], [No], [Poor ($L_2$)], [Semi-supervised GPS / IR],
    [*21*], [MacQueen's K-Means (1967)], [$O(n k d T)$], [$O(k d)$ online], [Order], [Sequential $L_2$ Voronoi], [Order-sensitive], [No], [Poor ($L_2$)], [Online single-pass stream],
    [*22*], [Forgy's K-Means (1965)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Batch $L_2$ Voronoi], [Fragile ($1 / n$)], [No], [Poor ($L_2$)], [Classical batch baseline],
    [*23*], [K-Means|| (2012)], [$O(n d r ell + r ell k d)$], [$O(n d / P + r ell d)$], [Seed], [Parallel $D^2$ Voronoi], [Fragile ($1 / n$)], [No], [Fair (dist.)], [Spark / distributed cluster],
    [*24*], [FastPAM / FasterPAM (2019)], [$O((n - k)^2 T)$], [$O(n^2)$], [Det.], [Exact metric medoids], [Robust (medoids)], [No], [Fair (metric)], [High-$k$ exact medoids],
    [*25*], [CLARA-MED (1990 / Ext.)], [$O(q (s^2 T + n k))$], [$O(s^2)$ stream], [Seed], [Median-scored medoids], [Very robust (med.)], [No], [Fair (metric)], [Corrupted large metric DB],
    [*26*], [Trimmed K-Means (1997)], [$O(n (k d + log n) T)$], [$O((n + k) d)$], [Seed], [Impartial trimmed $L_2$], [Explicit $alpha$-trim], [No], [Poor ($L_2$)], [Contaminated tabular data],
    [*27*], [Power K-Means (2019)], [$O(n k d T)$], [$O(n (d + k))$], [Seed], [Annealed power-$L_2$], [Moderate (homotopy)], [No], [Good (anneal)], [High-$d$ non-convex init.],
    [*28*], [W-K-Means (2005)], [$O(n k d T)$], [$O((n + k) d)$], [Seed], [Feature-weighted $L_2$], [Downweights noisy $d$], [No], [Good (wt.)], [Noisy feature tabular data],
    [*29*], [Minkowski W-K-Means (2012)], [$O(n d (k + T_"in") T)$], [$O((n + k) d)$], [Seed], [Cluster-wise weighted $L_p$], [Robust ($p approx 1$)], [No], [Good (wt-$L_p$)], [Subspace-weighted tabular],
    [*30*], [LBG Algorithm (1980)], [$O(n k d T_"max")$], [$O((n + k) d)$], [Det.], [Dyadic codebook Voronoi], [Fragile ($1 / n$)], [No], [Fair (VQ)], [Speech / image VQ codecs],
  )
})

#subsec("B", "Table V — Cat. II: Density & Mode-Seeking (Algorithms 31–55)")

Category II comprises twenty-five nonparametric density-level-set and kernel mode-seeking engines. Every single algorithm (25 / 25, 100%) determines cluster count automatically without specifying $k$, and 24 / 25 provide explicit noise excision $C_(-1)$ or isolated density-mode separation (all except Blurring Mean Shift's over-blurring regime). Twelve spatial-index methods—led by DBSCAN (Alg. 31 @ester1996dbscan), HDBSCAN (Alg. 32 @campello2013hdbscan), OPTICS (Alg. 33 @ankerst1999optics), DENCLUE (Alg. 35 @hinneburg1998denclue), and DeBaCl (Alg. 46 @kent2013debacl)—run in $O(n log n)$ in low dimensions but degrade to $O(n^2)$ under high-dimensional range-query collapse @beyer1999nearest, @schubert2017dbscanrevisited. Only SNN (Alg. 36 @ertoz2003snn), RNN-DBSCAN (Alg. 45 @bryant2018rnndbscan), and SNN-DPC (Alg. 52 @liu2018snndpc) achieve strong high-dimensional survival via rank-based neighborhood overlaps, while #box[DBSCAN++] (Alg. 47 @jang2019dbscanpp) and Mean Shift++ (Alg. 55 @jang2021meanshiftpp) cut time complexity to linear scaling.

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table V: Complete 10-Dimension Theoretical Reference Matrix — Cat. II: Density & Mode-Seeking (Algs. 31–55; $N = 25$, Auto-$k$: 25/25, High-$d$ Good/Strong: 3/25, Strictly Det.: 16/25).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*31*], [DBSCAN (1996)], [$O(n log n)$ / $O(n^2)$], [$O(n d)$], [Order], [Fixed-$epsilon$ level sets], [Native noise $C_(-1)$], [Yes ($epsilon$)], [Poor ($epsilon$-ball)], [GIS, spatial anomaly det.],
    [*32*], [HDBSCAN (2013/2015)], [$O(n log n)$ / $O(n^2)$], [$O(n d)$], [Det.], [Variable-density EoM], [Native noise $C_(-1)$], [Yes (EoM)], [Fair (mreach)], [Multi-density exploration],
    [*33*], [OPTICS (1999)], [$O(n log n)$ / $O(n^2)$], [$O(n d)$], [Det.], [Multi-scale valleys], [Native noise $C_(-1)$], [Yes ($xi$)], [Poor ($epsilon$-ball)], [Multi-scale density plots],
    [*34*], [Mean Shift (1995/2002)], [$O(n^2 d T)$], [$O(n d)$], [Det.], [KDE Morse mode basins], [Isolated outlier modes], [Yes ($h$)], [Poor (KDE)], [Image seg., visual track.],
    [*35*], [DENCLUE (1998)], [$O(n log n + n d T)$], [$O(n d)$], [Det.], [Analytical KDE attractors], [Native threshold $xi$], [Yes ($h, xi$)], [Poor (grid)], [Large scientific KDE data],
    [*36*], [SNN Clustering (2003)], [$O(n^2 d + n k^2)$], [$O(n k)$], [Det.], [Rank-overlap graph], [Native noise $C_(-1)$], [Yes (SNN)], [Strong (rank)], [High-$d$ climate, text, bio.],
    [*37*], [PDBSCAN (1999)], [$O((n / P) log n)$], [$O(n d / P)$], [Order], [Exact parallel DBSCAN], [Native noise $C_(-1)$], [Yes ($epsilon$)], [Poor ($R^*$-tree)], [Distributed spatial DBs],
    [*38*], [DBCLASD (1998)], [$O(n log n)$ / $O(n^2)$], [$O(n d)$], [Order], [Poisson NN-dist. basins], [Native noise], [Yes (auto)], [Poor (NN-dist)], [Parameter-free spatial GIS],
    [*39*], [GDBSCAN (1998)], [$O(n T_"pred")$], [$O(n d)$], [Order], [Predicate-defined sets], [Native noise $C_(-1)$], [Yes (pred.)], [Fair (pred.)], [Polygons, weighted graphs],
    [*40*], [ST-DBSCAN (2007)], [$O(n log n + n W_"max" d)$], [$O(n d)$], [Order], [Spatiotemporal cylinders], [Native ST noise], [Yes ($epsilon_1, epsilon_2$)], [Poor (ST)], [Event streams, meteorology],
    [*41*], [Quickshift (2008)], [$O(n^2 d)$], [$O(n d)$], [Det.], [Directed KDE parent tree], [Cut via $tau$ threshold], [Yes ($tau$)], [Poor (KDE)], [Fast image superpixels],
    [*42*], [Medoidshift (2007)], [$O(n^2 (d + n))$], [$O(n^2)$], [Det.], [Discrete medoid modes], [Robust (medoid KDE)], [Yes ($h$)], [Fair (metric)], [Non-Euclidean video/mesh],
    [*43*], [DENCLUE 2.0 (2007)], [$O(n^2 d T_"avg")$], [$O(n d)$], [Det.], [EM-step KDE attractors], [Native threshold $xi$], [Yes ($h, xi$)], [Poor (KDE)], [Exact continuous KDE],
    [*44*], [K-DBSCAN (2010)], [$O(n k d T + sum n_r^2 d)$], [$O(n d)$], [Seed], [Local-$epsilon$ Voronoi cells], [Local noise $C_(-1)$], [Yes (local)], [Poor ($L_2$)], [Spatially varied density],
    [*45*], [RNN-DBSCAN (2018)], [$O(n^2 d)$], [$O(n k)$], [Det.], [Reverse-$k$NN manifolds], [Native noise $C_(-1)$], [Yes ($k$)], [Good (RkNN)], [Multi-density metric data],
    [*46*], [DeBaCl (2013)], [$O(n k log n)$], [$O(n k)$], [Det.], [Level-set cluster tree], [Pruned low-mass leaves], [Yes (tree)], [Poor ($k$NN)], [Statistical mode trees],
    [*47*], [DBSCAN++ (2019)], [$O(m n d)$], [$O(n d)$], [Seed], [Subsampled core level sets], [Native noise $C_(-1)$], [Yes ($epsilon$)], [Fair ($m << n$)], [Large-scale density data],
    [*48*], [GriDBSCAN (2008)], [$O(n d + sum |tilde(P)_g|^2 d)$], [$O(n d)$], [Order], [Grid-boxed exact DBSCAN], [Native noise $C_(-1)$], [Yes ($epsilon$)], [Poor (grid)], [Low-$d$ spatial speedup],
    [*49*], [MR-DBSCAN (2011)], [$O((n \/ P) log(n \/ P))$], [$O((n d) \/ P)$], [Det.], [Cost-balanced MapReduce], [Native noise $C_(-1)$], [Yes ($epsilon$)], [Poor (grid)], [Hadoop / Spark spatial],
    [*50*], [VDBSCAN (2007)], [$O(n^2 + L n)$], [$O(n^2)$], [Order], [Multi-$epsilon$ peeled sets], [Native noise $C_(-1)$], [Yes ($k$-dist)], [Poor ($L_2$)], [Stratified multi-density],
    [*51*], [DPC (2014)], [$O(n^2)$], [$O(n^2)$], [Det.], [Decision-graph peaks], [Halo border excision], [Yes (graph)], [Fair ($d_c$)], [Biophysics, vision, tabular],
    [*52*], [SNN-DPC (2018)], [$O(n^2 + n k^2)$], [$O(n^2)$], [Det.], [SNN-weighted peaks], [Halo + rank resilience], [Yes (graph)], [Strong (SNN)], [High-$d$ / manifold peaks],
    [*53*], [OPTICS-OF (1999)], [$O(n log n)$ / $O(n^2)$], [$O(n d)$], [Det.], [Reachability + LOF score], [Explicit local outliers], [Yes ($xi$)], [Poor ($epsilon$-ball)], [Joint clustering + anomaly],
    [*54*], [Blurring Mean Shift (1995)], [$O(T n^2 d)$], [$O(n^2 + n d)$], [Det.], [Contracting KDE modes], [Sensitive to over-blur], [Yes ($h$)], [Poor (KDE)], [Fast cubic-conv. vision],
    [*55*], [Mean Shift++ (2021)], [$O(T n d 3^d)$], [$O((n + M_0) d)$], [Det.], [Grid-smoothed KDE modes], [Isolated grid cells], [Yes ($h$)], [Poor ($3^d$)], [Real-time 2D/3D/5D vision],
  )
})

#subsec("C", "Table VI — Cat. III: Hierarchical & Tree-Based (Algorithms 56–80)")

Category III catalogs twenty-five agglomerative, divisive, and incremental tree-structured algorithms producing nested multi-resolution dendrograms. Classical Lance–Williams linkages (Algs. 56–62, 72, 74 @lance1967general, @mullner2011modern) are strictly deterministic and complete in $O(n^2)$ via nearest-neighbor chains whenever reducibility holds (SLINK @sibson1973slink, CLINK @defays1977clink, UPGMA @sokal1958statistical, Ward @ward1963hierarchical), whereas DIANA (Alg. 63), Minimax Linkage (Alg. 73 @bien2011minimax), and Neighbor-Joining (Alg. 75 @saitou1987neighbor) pay $O(n^3)$. To escape the $O(n^2)$ RAM wall at massive scale, BIRCH (Alg. 64 @zhang1996birch), BUBBLE (Alg. 71 @ganti1999bubble), and coreset-backed BICO (Alg. 78 @fichtenberger2013bico) achieve single-pass $O(n)$ compression with sublinear RAM, while PERCH (Alg. 76 @kobren2017perch), GRINCH (Alg. 77 @monath2019grinch), SCC (Alg. 79 @monath2021scc), and parameter-free FINCH (Alg. 80 @sarfraz2019finch) scale hierarchical trees to millions of high-dimensional embeddings.

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table VI: Complete 10-Dimension Theoretical Reference Matrix — Cat. III: Hierarchical & Tree-Based (Algs. 56–80; $N = 25$, Auto-$k$: 19/25, High-$d$ Good/Strong: 7/25, Strictly Det.: 19/25).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*56*], [Single Linkage / SLINK (1973)], [$O(n^2)$], [$O(n)$ stream], [Det.], [MST sub-dom. ultrametric], [Chaining via bridges], [Dendro.], [Poor (chain)], [Topology, elongated paths],
    [*57*], [Complete Linkage / CLINK (1977)], [$O(n^2)$], [$O(n)$ stream], [Det.], [Compact clique diameter], [Outliers stay isolated], [Dendro.], [Fair (diam.)], [Compact bio. / taxonomy],
    [*58*], [Average Linkage / UPGMA (1958)], [$O(n^2)$], [$O(n^2)$], [Det.], [Mean inter-cluster dist.], [Moderate robustness], [Dendro.], [Fair (avg.)], [Genomics, phylogenetics],
    [*59*], [Ward's Method (1963)], [$O(n^2)$], [$O(n^2)$ / $O(n d)$], [Det.], [Min-variance spherical $L_2$], [Fragile (var. inflate)], [Dendro.], [Fair ($L_2$)], [General hierarchical bench.],
    [*60*], [WPGMA (1958)], [$O(n^2)$], [$O(n^2)$], [Det.], [Branch-equal ultrametric], [Moderate robustness], [Dendro.], [Fair (avg.)], [Unbalanced phylogeny],
    [*61*], [Centroid Linkage / UPGMC (1958)], [$O(n^2 log n)$], [$O(n^2)$], [Det.], [Centroid dist. (inversions)], [Fragile ($1 / n$)], [Dendro.], [Poor (inv.)], [Historical Euclidean tree],
    [*62*], [Median Linkage / WPGMC (1967)], [$O(n^2 log n)$], [$O(n^2)$], [Det.], [Midpoint dist. (inversions)], [Moderate (unweighted)], [Dendro.], [Poor (inv.)], [Size-invariant phylogeny],
    [*63*], [DIANA (1990)], [$O(n^3)$], [$O(n^2)$], [Det.], [Top-down max-diameter], [Splits outliers early], [Dendro.], [Fair (metric)], [Macro-split diagnostics],
    [*64*], [BIRCH (1996)], [$O(n d B log_B m)$], [$O(m d)$ CF-tree], [Order], [Spherical CF-ball tree], [Optional sparse prune], [No (macro)], [Poor ($L_2$)], [Massive DB pre-clustering],
    [*65*], [CURE (1998)], [$O(n_s^2 log n_s + n k c d)$], [$O(n_s d)$], [Seed], [Shrunk multi-rep. shapes], [Prunes slow-growing], [No], [Fair (rep.)], [Large non-spherical DBs],
    [*66*], [CHAMELEON (1999)], [$O(n^2 d + m^2 log m)$], [$O(k n)$], [Det.], [Dynamic RI + RC shapes], [Noise isolated in $k$NN], [No (cut)], [Good ($k$NN)], [Complex 2D/spatial shapes],
    [*67*], [ROCK (2000)], [$O(n^2 log n + n m_m m_a)$], [$O(n^2)$], [Det.], [Jaccard link-count sets], [Zero-link noise filter], [No], [Good (links)], [Market basket, categorical],
    [*68*], [MONA (1990)], [$O(n d^3)$], [$O(n d)$], [Det.], [Axis-split binary rules], [Sensitive to bit flip], [Dendro.], [Fair (binary)], [Interpretable rule trees],
    [*69*], [Genie (2016)], [$O(n^2)$], [$O(n^2)$], [Det.], [Gini-balanced MST shapes], [Resists chaining/noise], [Dendro.], [Fair (Gini)], [Balanced non-convex data],
    [*70*], [HCA-DBSCAN (2016)], [$O(n^2)$], [$O(n^2)$], [Det.], [Density-core hierarchy], [Native core noise cut], [Dendro.], [Fair (core)], [Noisy multi-scale shapes],
    [*71*], [BUBBLE (1999)], [$O(n (B log_B n_L + s) + n_L^2)$], [$O(n_L s)$], [Order], [Metric clustroid CF-tree], [Sparse leaf pruning], [No (macro)], [Fair (metric)], [Large non-Euclidean DBs],
    [*72*], [Lance-Williams (1967)], [$O(n^2)$–$O(n^2 log n)$], [$O(n^2)$], [Det.], [Parametric space warping], [Varies by $(alpha, beta, gamma)$], [Dendro.], [Fair (gen.)], [Unified linkage framework],
    [*73*], [Minimax Linkage (2011)], [$O(n^3)$], [$O(n^2)$], [Det.], [Prototype minimax balls], [Robust (medoid-like)], [Dendro.], [Fair (proto.)], [Interpretable bio. / med.],
    [*74*], [Flexible Beta Linkage (1967)], [$O(n^2)$ / $O(n^2 log n)$], [$O(n^2)$], [Det.], [Tunable dilation ($beta < 0$)], [Resists chaining ($beta < 0$)], [Dendro.], [Fair (tun.)], [Ecology, morphometrics],
    [*75*], [Neighbor-Joining (1987)], [$O(n^3)$], [$O(n^2)$], [Det.], [Additive metric tree], [Consistent under noise], [Dendro.], [Fair (add.)], [Molecular phylogenetics],
    [*76*], [PERCH (2017)], [$O(N d log N)$], [$O(N d)$], [Order], [Bounding-box rotations], [Collapse-mode bounding], [Dendro.], [Good (A\*)], [Extreme entity resolution],
    [*77*], [GRINCH (2019)], [$O(N d log N)$], [$O(N d)$], [Order], [Grafted cosine/linkage], [Graft fixes order noise], [Dendro.], [Strong (emb.)], [Extreme NLP / entity trees],
    [*78*], [BICO (2013)], [$O(N d log S + S^2 d)$], [$O(S d)$], [Order], [$(1+epsilon)$ $k$-means coreset], [Fragile ($L_2$ coreset)], [No (macro)], [Good (stream)], [Big-data guaranteed VQ],
    [*79*], [SCC (2021)], [$O(N d log N + R |E|)$], [$O(N d + |E|)$], [Det.], [Multi-round reciprocal NN], [Threshold isolates noise], [Dendro.], [Strong (ANN)], [Billion-scale embeddings],
    [*80*], [FINCH (2019)], [$O(n log n)$ (ANN)], [$O(n d)$], [Det.], [1-NN graph contraction], [Merges outliers late], [Yes (auto)], [Strong (1-NN)], [Parameter-free vision/emb.],
  )
})

#subsec("D", "Table VII — Cat. IV: Distribution, Statistical & Topic Models (Algorithms 81–102)")

Category IV unifies twenty-two generative finite mixtures, nonparametric Bayesian processes, and statistical topic factorizations. Twenty of twenty-two methods fall in the bilinear iterative EM / CAVI / collapsed-Gibbs tier ($O(T n k d)$ to $O(T n k d^2 + T k d^3)$), whereas incremental conceptual trees COBWEB (Alg. 85 @fisher1987cobweb) and CLASSIT (Alg. 86 @gennari1989classit) operate in linear single-pass time. Full-covariance continuous mixtures—GMM (Alg. 81 @dempster1977em), TMM (Alg. 88 @peel2000robust), and Skew-Normal (Alg. 98 @lin2009skewnormal)—incur $O(k d^2)$ parameter explosion when $d approx n$; parsimonious eigen-decompositions in MCLUST / Mixmod (Algs. 93–94 @fraley2002mclust, @biernacki2006mixmod), directional vMF mixtures (Alg. 89 @banerjee2005vmf), and discrete simplex models (LDA @blei2003lda, pLSA @hofmann1999plsa, HDP @teh2006hdp, NMF @lee1999nmf, DMM @yin2014dmm, CTM @blei2007ctm, BTM @yan2013btm) restore strong high-dimensional performance, while Student-$t$ mixtures (Alg. 88) provide bounded-influence outlier protection.

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table VII: Complete 10-Dimension Theoretical Reference Matrix — Cat. IV: Distribution, Statistical & Topic Models (Algs. 81–102; $N = 22$, Auto-$k$: 10/22, High-$d$ Good/Strong: 13/22, Strictly Det.: 1/22).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*81*], [Gaussian Mixture / GMM (1977)], [$O(T n k d^2 + T k d^3)$], [$O(n k + k d^2)$], [Seed], [Full ellipsoidal Gaussians], [Fragile ($1 / n$)], [No], [Poor ($d^2$)], [Density est., speech/audio],
    [*82*], [Bayesian GMM / DPGMM (2000)], [$O(T n k d^2 + T k d^3)$], [$O(n k + k d^2)$], [Seed], [NIW-regularized ellipsoids], [Prior prevents collapse], [Yes (DP)], [Fair (prior)], [Nonparametric model sel.],
    [*83*], [LDA (2003)], [$O(T_"out" (T_"in" k "nnz" + k V))$], [$O(k (V + n))$], [Stoch.], [Dirichlet topic simplex], [Smoothed by $alpha, beta$], [No], [Strong (BoW)], [Text corpora, topic mining],
    [*84*], [pLSA (1999)], [$O(T k "nnz"(X))$], [$O(k (n + V))$], [Seed], [Multinomial KL simplex], [Overfits unseen docs], [No], [Strong (BoW)], [IR, semantic indexing],
    [*85*], [COBWEB (1987)], [$O(n B^2 d V h)$], [$O(|cal(T)| d V)$], [Order], [Category-utility concepts], [Split/merge absorbs], [Dendro.], [Fair (cat.)], [Cognitive science, symbolic],
    [*86*], [CLASSIT (1989)], [$O(n B^2 d h)$], [$O(|cal(T)| d)$], [Order], [Continuous Gaussian CU], [Acuity floor $sigma_"min"$], [Dendro.], [Fair (diag.)], [Incremental numerical tax.],
    [*87*], [AutoClass (1996)], [$O(R T n k (d_c + d_m V))$], [$O(n k + k d V)$], [Seed], [Mixed indep. distributions], [Moderate (Bayes MAP)], [Yes (Bayes)], [Good (indep.)], [Astronomy, NASA mixed DB],
    [*88*], [Student's t-Mixture / TMM (2000)], [$O(T n k d^2 + T k d^3)$], [$O(n k + k d^2)$], [Seed], [Heavy-tailed ellipsoids], [Very robust (bounded IF)], [No], [Poor ($d^2$)], [Finance, medical outliers],
    [*89*], [vMF Mixture (2005)], [$O(T n k d)$], [$O((n + k) d + n k)$], [Seed], [Directional caps on $bb(S)^(d-1)$], [Moderate (unit sphere)], [No], [Strong (cos)], [Text, NLP, gene expression],
    [*90*], [Latent Class Analysis (1968)], [$O(T k (n d + M))$], [$O(n d + n k + k M)$], [Seed], [Cond.-indep. categorical], [Moderate (discrete)], [No], [Good (cat.)], [Psychometrics, epidemiology],
    [*91*], [HDP (2006)], [$O(T n k)$], [$O(n + J k + k d)$], [Stoch.], [Shared grouped topics], [Hierarchical smoothing], [Yes (HDP)], [Strong (BoW)], [Multi-corpus topic discovery],
    [*92*], [NMF (1999)], [$O(T n d k)$], [$O("nnz"(X) + (n + d) k)$], [Seed], [Additive parts-based cone], [KL/L1 variants robust], [No], [Strong (low-$k$)], [Spectra, vision, genomics],
    [*93*], [MCLUST (2002)], [$O(M T k (n d^2 + d^3))$], [$O(n k + k d^2)$], [Det.], [14 eigen-decomposed GMMs], [Optional noise comp.], [Yes (BIC)], [Fair (pars.)], [Biostatistics, R gold std.],
    [*94*], [Mixmod (2006)], [$O(M T k (n d^2 + d^3))$], [$O(n k + k d^2)$], [Seed], [28 Gauss./multinom. models], [Parsimonious reg.], [Yes (BIC/ICL)], [Fair (pars.)], [Mixed tabular / clinical],
    [*95*], [Poisson Mixture (1977)], [$O(T n k d)$], [$O(n k + k d)$], [Seed], [I-divergence count cone], [Overdispersion-sens.], [No], [Good (counts)], [RNA-seq, click/event counts],
    [*96*], [Bernoulli Mixture (1968)], [$O(T n k d)$], [$O(n k + k d)$], [Seed], [Binary hypercube prod.], [Moderate (bounded)], [No], [Good (bin.)], [Binary docs, digits, logs],
    [*97*], [Dirichlet Multinom. / DMM (2014)], [$O(T k bar(L) n)$], [$O("nnz"(X) + k V)$], [Stoch.], [Single-topic short text], [Symmetric Dir. prior], [Yes (collapse)], [Strong (text)], [Tweets, search queries],
    [*98*], [Skew-Normal Mixture (2009)], [$O(T k (n d^2 + d^3))$], [$O(n k + k d^2)$], [Seed], [Asymmetric skewed hulls], [Sensitive (light tail)], [No], [Poor ($d^2$)], [Flow cytometry, income],
    [*99*], [Infinite GMM / IGMM (2000)], [$O(T k_"eff" n d^2)$], [$O(n d + k_"eff" d^2)$], [Stoch.], [CRP nonparametric GMM], [Spawns outlier clusters], [Yes (CRP)], [Fair (Bayes)], [Unknown-$k$ Bayes density],
    [*100*], [Correlated Topic / CTM (2007)], [$O(T n (k bar(L) + k^2) + T k^3)$], [$O(k V + k^2 + n k)$], [Seed], [Logistic-normal simplex], [Cov. captures co-occ.], [No], [Strong (BoW)], [Correlated scientific text],
    [*101*], [Biterm Topic / BTM (2013)], [$O(T k |cal(B)|)$], [$O(|cal(B)| + k V)$], [Stoch.], [Global word-pair simplex], [Robust to short sparsity], [No], [Strong (text)], [Short social posts, ads],
    [*102*], [SNOB (1968)], [$O(T n k d)$], [$O(n d + n k + k d)$], [Seed], [MML diagonal/mixed mix.], [MML prunes spurious], [Yes (MML)], [Good (diag.)], [MDL / MML statistical DBs],
  )
})

#subsec("E", "Table VIII — Cat. V: Graph, Spectral & Network (Algorithms 103–132)")

Category V spans thirty spectral cut relaxations, modularity/flow community detectors, multilevel partitioners, overlapping cover extractors, and stochastic block models. Its complexity profile is starkly bimodal: dense spectral eigensolvers (NJW @ng2002spectral, Ncut @shi2000ncut, Ratio Cut @hagen1992ratiocut), Girvan–Newman (Alg. 112 @girvan2002community), HCS (Alg. 118 @hartuv2000hcs), and Kluger Spectral Biclustering (Alg. 123 @kluger2003spectral) require $O(n^3)$ or $O(m^2 n)$ exact computation, whereas sparse CSR message-passing and multilevel refinement—Louvain (Alg. 107 @blondel2008louvain), Leiden (Alg. 108 @traag2019leiden), LPA (Alg. 110 @raghavan2007lpa), Chinese Whispers (Alg. 111 @biemann2006chinese), Infomap (Alg. 113 @rosvall2008infomap), FluidC (Alg. 119 @pares2018fluidc), METIS (Alg. 121 @karypis1998metis), SLPA (Alg. 126 @xie2011slpa), and PIC (Alg. 132 @lin2010pic)—achieve near-linear $O(m)$ or $O(m log n)$ time and $O(n + m)$ RAM on billion-edge graphs. Crucially, Leiden fixes Louvain's disconnected-community defect, DC-SBM (Alg. 131 @karrer2011dcsbm) resolves scale-free hub bias, and BigCLAM / OSLOM / CPM / DEMON (Algs. 124–128) enable overlapping communities.

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table VIII: Complete 10-Dimension Theoretical Reference Matrix — Cat. V: Graph, Spectral & Network (Algs. 103–132; $N = 30$, Auto-$k$: 19/30, High-$d$ Good/Strong: 29/30, Strictly Det.: 12/30).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*103*], [Spectral NJW (2002)], [$O(n^2 d + n^3)$], [$O(n^2)$], [Seed], [Normalized $L_"sym"$ sphere], [Robust (deg-norm.)], [No (eigengap)], [Good ($k$NN)], [Non-convex manifolds, graphs],
    [*104*], [Normalized Cuts / Ncut (2000)], [$O(n^3)$ (exact)], [$O(n^2)$], [Det.], [Volume-balanced cuts], [Resists isolated nodes], [No], [Good (graph)], [Image segmentation, vision],
    [*105*], [Ratio Cut (1992)], [$O(n^3)$ (exact)], [$O(n^2)$], [Det.], [Cardinality-balanced $L$], [Sensitive to low-deg.], [No], [Good (graph)], [VLSI / circuit partitioning],
    [*106*], [Affinity Propagation (2007)], [$O(T n^2)$], [$O(n^2)$], [Det.], [Exemplar factor graph], [Exemplar-robust], [Yes ($s(i,i)$)], [Fair (sim.)], [Face images, exon arrays],
    [*107*], [Louvain Method (2008)], [$O(L m)$], [$O(n + m)$], [Order], [Multilevel modularity $Q$], [Res. limit + disconn.], [Yes ($Q$)], [Strong (sparse)], [Web & social networks],
    [*108*], [Leiden Algorithm (2019)], [$O(L m)$], [$O(n + m)$], [Stoch.], [Connected CPM / $Q$], [Guaranteed connected], [Yes ($Q$/CPM)], [Strong (sparse)], [scRNA-seq, huge networks],
    [*109*], [Markov Clustering / MCL (2000)], [$O(T n k^2)$ (pruned)], [$O(n k)$], [Det.], [Random-walk flow traps], [Isolated nodes pruned], [Yes (infl. $r$)], [Good (sparse)], [Protein interaction (PPI)],
    [*110*], [Label Propagation / LPA (2007)], [$O(T m)$], [$O(n + m)$], [Stoch.], [Local majority-vote], [Monster-cluster risk], [Yes (emerg.)], [Strong (sparse)], [Real-time huge graphs],
    [*111*], [Chinese Whispers (2006)], [$O(T m)$], [$O(n + m)$], [Stoch.], [Weighted local majority], [Tolerates edge noise], [Yes (emerg.)], [Strong (sparse)], [NLP word-sense graphs],
    [*112*], [Girvan-Newman (2002)], [$O(m^2 n)$], [$O(n + m)$], [Det.], [Edge-betweenness bridges], [Noise bridges cut first], [Dendro. ($Q$)], [Good (small)], [Small social/citation nets],
    [*113*], [Infomap (2008)], [$O(m log n)$], [$O(n + m)$], [Stoch.], [Directed flow MDL], [Flow-robust (no res. lim.)], [Yes (MDL)], [Strong (flow)], [Citation, web, transport],
    [*114*], [Walktrap (2006)], [$O(m n H)$], [$O(n^2)$], [Det.], [$t$-step diffusion Ward], [Short walks smooth noise], [Dendro. ($Q$)], [Good (graph)], [Medium complex networks],
    [*115*], [Fast Greedy / CNM (2004)], [$O(m d log n)$], [$O(n + m)$], [Det.], [Agglom. modularity $Q$], [Skewed hub merging], [Dendro. ($Q$)], [Good (sparse)], [Large sparse social nets],
    [*116*], [Leading Eigenvector (2006)], [$O(n (m + n))$], [$O(n + m)$], [Det.], [Spectral modularity $B$], [Moderate (global $B$)], [Yes ($Delta Q$)], [Good (sparse)], [Spectral community splits],
    [*117*], [Spinglass (2006)], [$O(N_"sw" (m + n q))$], [$O(n + m + q)$], [Stoch.], [Potts Hamiltonian ($gamma$)], [Supports signed edges], [Yes ($gamma$)], [Good (signed)], [Signed / physics networks],
    [*118*], [HCS (2000)], [$O(N_"cut" (m n + n^2 log n))$], [$O(n + m)$], [Det.], [$> n / 2$ edge connectivity], [Drops sparse singletons], [Yes (min-cut)], [Good (dense)], [Gene expression / cDNA],
    [*119*], [FluidC (2018)], [$O(T (n + m))$], [$O(n + m)$], [Stoch.], [Density-inhibited fluids], [Requires connected $G$], [No], [Strong (sparse)], [Fixed-$k$ scalable graphs],
    [*120*], [Kernighan-Lin (1970)], [$O(P n^2 log n)$], [$O(n + m)$], [Order], [Balanced 2-way min-cut], [Strict balance constraint], [No ($k=2$)], [Good (graph)], [VLSI / parallel mesh split],
    [*121*], [METIS (1998)], [$O((n + m) log k)$], [$O(n + m)$], [Order], [Multilevel balanced $k$-cut], [Coarsening smooths noise], [No], [Strong (CSR)], [HPC FEM mesh / sparse mat.],
    [*122*], [Chaco (1995)], [$O((m + n log n) log k)$], [$O(n + m)$], [Det.], [Multilevel spectral cut], [Lanczos smooths noise], [No], [Strong (CSR)], [Scientific PDE domain dec.],
    [*123*], [Spectral Biclustering (2003)], [$O(n m min(n, m))$], [$O(n m)$], [Det.], [Checkerboard bistoch. SVD], [Sensitive to raw scale], [No], [Strong (SVD)], [Microarray / gene matrices],
    [*124*], [OSLOM (2011)], [$O(n^2)$], [$O(n + m)$], [Stoch.], [Order-stat. $p$-value cover], [Explicit homeless nodes], [Yes ($p$-val)], [Good (stat.)], [Overlapping noisy networks],
    [*125*], [BigCLAM (2013)], [$O(T k (m + n))$], [$O(n k + m)$], [Seed], [Dense-overlap NMF affil.], [Background prob. $epsilon$], [Yes (holdout)], [Strong (sparse)], [Million-node overlapping],
    [*126*], [SLPA (2011)], [$O(T m)$], [$O(T n + m)$], [Stoch.], [Memory-bag overlap LPA], [Threshold $r$ drops noise], [Yes (emerg.)], [Strong (sparse)], [Overlapping social nets],
    [*127*], [CPM (2005)], [$O(d n 3^(d / 3) + c_k^2 k)$], [$O(n + m + c_k)$], [Det.], [Rolling $k$-clique cover], [Drops non-clique nodes], [Yes ($k$-clq)], [Good (sparse)], [Protein complexes, social],
    [*128*], [DEMON (2012)], [$O(T tau(G) + R M bar(s))$], [$O(n + m)$], [Stoch.], [Ego-network local cover], [Local ego filters noise], [Yes (ego)], [Strong (ego)], [Ego-centric social graphs],
    [*129*], [Fast Unfolding (2008)], [$O(L T m)$], [$O(n + m)$], [Order], [Exact Louvain formaliz.], [Resolution-limit bound], [Yes ($Q$)], [Strong (sparse)], [Telco / web billion-edge],
    [*130*], [SBM (1983)], [$O(T (m k + n k^2))$], [$O(n k + m)$], [Seed], [Assort./disassort. blocks], [Principled edge likelihood], [No (MDL ext.)], [Strong (gen.)], [Network science inference],
    [*131*], [Degree-Corrected SBM (2011)], [$O(T k (m + n))$], [$O(n k + m)$], [Seed], [Hub-tolerant block graph], [Decouples hubs vs. blocks], [No (MDL ext.)], [Strong (hubs)], [Scale-free real networks],
    [*132*], [Power Iteration / PIC (2010)], [$O(T m + I n k)$], [$O(n + m)$], [Seed], [1D deflated RW subspace], [Moderate (power smooth)], [No], [Good (sparse)], [Fast large-graph spectral],
  )
})

#subsec("F", "Table IX — Cat. VI: Grid, Subspace & Projection (Algorithms 133–150)")

Category VI contains eighteen grid-quantized and high-dimensional subspace clustering methods explicitly engineered to defeat Beyer's curse of dimensionality @beyer1999nearest, @kriegel2009subspacesurvey—achieving the highest continuous-feature high-$d$ resilience share of any Euclidean family (15 / 18 rated Good or Strong). Low-dimensional grid engines STING (Alg. 133 @wang1997sting) and WaveCluster (Alg. 135 @sheikholeslami1998wavecluster) run in $O(n + K)$ linear time with $O(K)$ memory but collapse exponentially in $d$. Conversely, exhaustive bottom-up Apriori subspace finders—CLIQUE (Alg. 134 @agrawal1998clique), SUBCLU (Alg. 137 @kailing2004subclu), DOC (Alg. 139 @procopiuc2002doc), MAFIA (Alg. 140 @goil1999mafia), and ENCLUS (Alg. 143 @cheng1999enclus)—pay worst-case exponential $O(2^(k_"max"))$ lattice search. Polynomial-time projected and statistical alternatives—PROCLUS (Alg. 136 @aggarwal1999proclus), ORCLUS (Alg. 138 @aggarwal2000orclus), OptiGrid (Alg. 141 @hinneburg1999optigrid), P3C (Alg. 148 @moise2006p3c), and STATPC (Alg. 149 @moise2008statpc)—avoid the exponential barrier, while Sparse Subspace Clustering (Alg. 150 @elhamifar2013ssc) recovers arbitrarily oriented oblique subspaces with $ell_1$ sparse-error guarantees at $O(n^3)$ spectral cost.

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table IX: Complete 10-Dimension Theoretical Reference Matrix — Cat. VI: Grid, Subspace & Projection (Algs. 133–150; $N = 18$, Auto-$k$: 14/18, High-$d$ Good/Strong: 15/18, Strictly Det.: 11/18).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*133*], [STING (1997)], [$O(n + K)$], [$O(K)$ cells], [Det.], [Hierarchical rect. cells], [Cell stats filter outliers], [Yes (query)], [Poor ($c^d$)], [Spatial GIS SQL queries],
    [*134*], [CLIQUE (1998)], [$O(c^(k_"s") + m k_"s" n)$], [$O(n d + |cal(U)|)$], [Det.], [Axis-parallel dense boxes], [Native low-density cut], [Yes (grid)], [Good (low-$s$)], [High-$d$ axis subspaces],
    [*135*], [WaveCluster (1998)], [$O(n + K)$], [$O(K)$ grid], [Det.], [Multi-res. wavelet shapes], [Wavelet low-pass denoise], [Yes (scale)], [Poor ($c^d$)], [Large 2D/3D spatial DBs],
    [*136*], [PROCLUS (1999)], [$O(I n k d + k d log(k d))$], [$O(n d)$], [Seed], [Disjoint axis medoids], [Explicit outlier set], [No], [Strong (proj.)], [High-$d$ tabular / finance],
    [*137*], [SUBCLU (2004)], [$O(|cal(S)_"cand"| n^2)$], [$O(n^2)$], [Det.], [Arbitrary-shape subspace], [Native DBSCAN noise], [Yes ($epsilon$)], [Good (subsp.)], [Non-convex subspace data],
    [*138*], [ORCLUS (2000)], [$O(k_0^3 + k_0 n d + k_0^2 d^3)$], [$O(n d + k_0 d^2)$], [Seed], [Oblique PCA subspaces], [Outlier distance prune], [No], [Strong (PCA)], [Correlated high-$d$ features],
    [*139*], [DOC (2002)], [$O(n d alpha^(-(r+1)))$], [$O(n d)$], [Stoch.], [Monte Carlo hyper-boxes], [Explicit unclustered set], [Yes (iter.)], [Strong (MC)], [Theoretical proj. clusters],
    [*140*], [MAFIA (1999)], [$O(n (d + |cal(U)_"cand"|))$], [$O(n d + |cal(U)|)$], [Det.], [Adaptive variable boxes], [Adaptive density threshold], [Yes (grid)], [Good (adapt.)], [Parallel scientific tables],
    [*141*], [OptiGrid (1999)], [$O(n d log n)$], [$O(n d)$], [Det.], [Optimal valley cut planes], [Sparse cells discarded], [Yes (cuts)], [Good (proj.)], [High-$d$ KDE valley search],
    [*142*], [BANG-Clustering (1998)], [$O(n d log n)$], [$O(n d)$], [Order], [Multi-scale BANG-tree], [Low-density tail drop], [Yes (dendro)], [Fair (grid)], [Spatial database indexing],
    [*143*], [ENCLUS (1999)], [$O(n d dot.op 2^(k_"max"))$], [$O(n d)$], [Det.], [Low-entropy subspaces], [Entropy screens uniform], [Yes (ent.)], [Good (ent.)], [Entropy feature discovery],
    [*144*], [FINDIT (2004)], [$O(n k d + |S| |M| d)$], [$O(n d)$], [Seed], [DOD-metric subspaces], [Explicit outlier class], [No], [Strong (DOD)], [High-$d$ attribute subsets],
    [*145*], [COSA (2004)], [$O(T n^2 d)$], [$O(n^2 + n d)$], [Det.], [Pairwise soft-subspace], [Outlier-resistant weights], [Dendro.], [Strong (wt.)], [Genomics, microarray profiles],
    [*146*], [FIRES (2005)], [$O(d n log n + m n^2 d_"avg")$], [$O(n d + m n)$], [Det.], [Merged 1D base clusters], [Base cluster noise filter], [Yes (merge)], [Strong (1D)], [Fast high-$d$ subspace DB],
    [*147*], [INSCY (2008)], [$O(n d log n + |cal(T)_"vis"| n^2)$], [$O(n d + n^2)$], [Det.], [Non-redundant dense sub.], [Density + redundancy cut], [Yes (idx)], [Strong (non-red)], [Redundancy-free subspaces],
    [*148*], [P3C (2006)], [$O(d n + |cal(K)| n + T n k d)$], [$O(n d)$], [Seed], [Poisson-core EM hulls], [Poisson $p$-val outlier set], [Yes (Pois.)], [Strong (1D-EM)], [Statistical high-$d$ data],
    [*149*], [STATPC (2008)], [$O(n_s d^2 n + |cal(H)|^2 n)$], [$O(n d)$], [Stoch.], [Non-red. stat. hyper-rect], [Explicit stat. outlier set], [Yes (stat.)], [Strong (stat.)], [High-$d$ non-redundant DB],
    [*150*], [Sparse Subspace / SSC (2013)], [$O(T min(n, d) n^2 + n^3)$], [$O(n^2)$], [Det.], [Union of oblique subspaces], [$ell_1$ sparse error + $L_2$], [No (eigengap)], [Strong ($ell_1$)], [Motion seg., face lighting],
  )
})

#subsec("G", "Table X — Cat. VII: Fuzzy, Biclustering, Time-Series & Domain (Algorithms 151–175)")

Category VII brings together twenty-five specialized engines across five applied sub-families: fuzzy/possibilistic soft partitioning (Algs. 151–156 @bezdek1981fcm, @krishnapuram1993pcm), two-mode gene-expression biclustering (Algs. 157–161 @cheng2000biclustering, @hochreiter2010fabia), elastic/shape time-series clustering (Algs. 162–165 @paparrizos2015kshape, @petitjean2011dba, @cuturi2017softdtw, @dempster2020rocket), computer-vision/LiDAR/trajectory engines (Algs. 166–169 @achanta2012slic, @vincent1991watersheds, @rusu2011pcl, @lee2007traclus), and ensemble/multi-view consensus methods (Algs. 170–175 @monti2003consensus, @strehl2002cluster, @cai2013multiview, @kumar2011coregularized). Standard FCM (Alg. 151) suffers the sum-to-one noise trap and high-$d$ fuzzification collapse ($u_(i j) -> 1 / k$), cured respectively by PCM (Alg. 152) and KFCM (Alg. 155 @zhang2004kfcm). Among time-series engines, k-Shape and ROCKET-Clustering avoid the $O(n m^2)$ quadratic warp bottleneck of DBA and Soft-DTW, while HGPA and MCLA (Algs. 172–173) eliminate the $O(n^2)$ co-association RAM wall of Monti Consensus and CSPA (Algs. 170–171).

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table X: Complete 10-Dimension Theoretical Reference Matrix — Cat. VII: Fuzzy, Biclustering, Time-Series & Domain (Algs. 151–175; $N = 25$, Auto-$k$: 7/25, High-$d$ Good/Strong: 14/25, Strictly Det.: 5/25).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*151*], [Fuzzy C-Means / FCM (1981)], [$O(n k d T)$], [$O(n (d + k))$], [Seed], [Soft $L_2$ power simplex], [Fragile (sum-to-1 trap)], [No], [Poor ($u -> 1 / k$)], [Medical MRI, soft control],
    [*152*], [Possibilistic / PCM (1993)], [$O(n k d T)$], [$O(n (d + k))$], [Seed], [Decoupled typicality], [Very robust (typicality)], [No], [Poor ($L_2$)], [Noisy fuzzy / sensor data],
    [*153*], [Gustafson-Kessel / GK (1979)], [$O((n k d^2 + k d^3) T)$], [$O(n k + k d^2)$], [Seed], [Fixed-vol. fuzzy ellipsoids], [Sensitive (cov. pull)], [No], [Poor ($d^2$)], [Elongated fuzzy control],
    [*154*], [Gath-Geva / GG (1989)], [$O(T k (n d^2 + d^3))$], [$O(n k + k d^2)$], [Seed], [MLE exp.-dist. ellipsoids], [Sensitive (exp. blowup)], [No], [Poor ($d^2$)], [Variable-density fuzzy],
    [*155*], [Kernel FCM / KFCM (2004)], [$O(T k n^2)$], [$O(n^2)$], [Seed], [RKHS non-convex fuzzy], [RBF kernel robust], [No], [Good (RKHS)], [Non-convex fuzzy / MRI],
    [*156*], [Spatial FCM / sFCM (2006)], [$O(T k n (d + |cal(N)|))$], [$O(n (d + k))$], [Seed], [Spatially smoothed fuzzy], [Suppresses salt-pepper], [No], [Poor (spatial)], [Brain MRI / remote imaging],
    [*157*], [Cheng-Church (2000)], [$O(K n p (n + p))$], [$O(n p)$], [Stoch.], [Low-MSR additive submat.], [Sensitive to mask noise], [No ($K$)], [Strong (2-mode)], [Gene expression microarray],
    [*158*], [PLAID (2002)], [$O(K T n p)$], [$O(n p)$], [Seed], [Overlapping additive layers], [Perm.-test noise screen], [Yes (perm.)], [Strong (2-mode)], [Transcriptomics, yeast],
    [*159*], [Bimax (2006)], [$O(n m beta min(n, m))$], [$O(n m)$], [Det.], [Exact all-1s submatrices], [Fragile (single 0 breaks)], [Yes (exact)], [Good (binary)], [Binary gene / regulatory],
    [*160*], [ISA (2003)], [$O(S T n m)$], [$O(n m)$], [Seed], [Thresholded co-reg. modules], [Z-score threshold robust], [Yes (seeds)], [Strong (2-mode)], [Co-regulated gene modules],
    [*161*], [FABIA (2010)], [$O(T (p^3 (n + m) + p n m))$], [$O(n m + p(n + m))$], [Seed], [Sparse multiplicative FA], [Heavy-tail Laplace robust], [No ($p$)], [Strong (sparse)], [Drug discovery, transcript.],
    [*162*], [k-Shape (2015)], [$O(T (n k m log m + k m^3))$], [$O(n m + m^2)$], [Seed], [Shift/scale-inv. SBD shape], [Moderate (z-normed)], [No], [Good (FFT)], [ECG, IoT, UCR time series],
    [*163*], [DBA (2011)], [$O(T (k + I) n m^2)$], [$O(n m + m^2)$], [Seed], [Elastic DTW barycenter], [Sensitive (warped mean)], [No], [Fair (band)], [Speech, gesture, ECG warp],
    [*164*], [Soft-DTW (2017)], [$O(T (k + I) n m^2)$], [$O(n m + m^2)$], [Seed], [Smooth $C^infinity$ warp], [Smoothed by $gamma > 0$], [No], [Fair (diff.)], [Deep / differentiable TS],
    [*165*], [ROCKET-Clustering (2020)], [$O(n M m + T k n M)$], [$O(n M)$], [Seed], [Random conv. feature $L_2$], [PPV features noise-res.], [No], [Strong (conv.)], [Fast long time-series],
    [*166*], [SLIC (2012)], [$O(N I)$], [$O(N)$], [Det.], [Compact 5D CIELAB+XY], [Local $2S times 2S$ cap], [No ($K$ spx)], [Poor (5D)], [Real-time superpixels],
    [*167*], [Watershed (1991)], [$O(N log N)$], [$O(N)$], [Det.], [Topographic basins], [Over-seg. unless markers], [Yes (minima)], [Poor (2D/3D)], [Cell/microscopy, CT seg.],
    [*168*], [DBSCAN-3D / PCL (2011)], [$O(N log N)$], [$O(N)$], [Det.], [3D Euclidean components], [Min/max-size outlier cut], [Yes ($epsilon$)], [Poor (3D $k$d)], [Autonomous LiDAR, robotics],
    [*169*], [TRACLUS (2007)], [$O(M^2)$], [$O(N_"pts" + M)$], [Order], [MDL sub-trajectory bundles], [Native segment noise], [Yes ($epsilon$)], [Poor (2D/3D)], [Hurricanes, animal/GPS],
    [*170*], [Consensus Clustering (2003)], [$O(H (T_"base" + n^2))$], [$O(n^2)$], [Stoch.], [Inherited co-assoc. tree], [Bootstrap variance smooth], [Yes (CDF)], [Good (meta)], [Cancer subtyping, omics],
    [*171*], [CSPA (2002)], [$O(n^2 (r + k))$], [$O(n^2)$], [Det.], [Co-association graph cut], [Ensemble consensus smooth], [No], [Good (ens.)], [Small/med. cluster ensemb.],
    [*172*], [HGPA (2002)], [$O(n k r)$], [$O(n r)$], [Order], [Balanced hypergraph cut], [Sensitive to imbalance], [No], [Good (ens.)], [Large balanced ensembles],
    [*173*], [MCLA (2002)], [$O(n k^2 r^2)$], [$O(n k + k^2 r^2)$], [Order], [Meta-cluster Jaccard graph], [Resists outlier runs], [No], [Good (ens.)], [Scalable soft ensembles],
    [*174*], [Multi-View K-Means (2013)], [$O(T n k D)$], [$O(n D + k D)$], [Seed], [Weighted multi-view $L_2$], [$ell_(2,1)$ view-noise rob.], [No], [Strong (views)], [Multimodal text+image],
    [*175*], [Co-reg. Multi-View Spec. (2011)], [$O(T V n^3)$], [$O(V n^2)$], [Seed], [Consistent view eigenspace], [Cross-view regularization], [No], [Strong (views)], [Multi-omics, multimedia],
  )
})

#subsec("H", "Table XI — Cat. VIII: Neural, Deep, Streaming & Evolutionary (Algorithms 176–200)")

Category VIII integrates twenty-five modern computational paradigms: competitive topological neural maps (Algs. 176–180 @kohonen1982som, @fritzke1995gng, @carpenter1987art), deep representation-clustering networks (Algs. 181–192 @xie2016dec, @jiang2017vade, @ji2019iic, @vangansbeke2020scan, @caron2020swav, @li2021contrastive), bounded-memory single-pass data-stream engines (Algs. 193–197 @aggarwal2003clustream, @cao2006denstream, @ackermann2012streamkmpp, @kranen2011clustree, @chen2007dstream), and stochastic population metaheuristics (Algs. 198–200 @krishna1999gka, @vandermerwe2003pso, @shelokar2004aco). Zero methods in Category VIII are strictly order-invariant deterministic: 16 / 25 are intrinsically stochastic (mini-batch SGD or evolutionary Markov search) and 7 / 25 are order-dependent online streams. Remarkably, Category VIII spans both extremes of the space-vs.-expressivity trade-off: ten neural/streaming engines run in strictly sublinear $o(n)$ RAM (e.g., CluStream, DenStream, StreamKM++, GNG, ART), while twelve deep perceptual networks achieve Strong high-dimensional recovery on raw vision, audio, and language via learned nonlinear encoders $Phi_omega$ regularized by reconstruction (IDEC, DCN), VAE-GMM priors (VaDE), mutual information (IIC, CC), or Sinkhorn optimal transport (SwAV).

#block(breakable: true, above: 0.45em, below: 0.85em, {
  set par(first-line-indent: 0pt, justify: false)
  text(size: 7.8pt, weight: "bold", fill: ieee_navy)[Table XI: Complete 10-Dimension Theoretical Reference Matrix — Cat. VIII: Neural, Deep, Streaming & Evolutionary (Algs. 176–200; $N = 25$, Auto-$k$: 5/25, High-$d$ Good/Strong: 15/25, Strictly Det.: 0/25).]
  v(3pt)
  set text(size: 6.75pt)
  show math.equation: set text(size: 6.75pt)
  table(
    columns: (19pt, 1.24fr, 1.15fr, 0.86fr, 0.37fr, 1.12fr, 0.85fr, 0.54fr, 0.56fr, 1.02fr),
    align: (center, left, left, left, center, left, left, center, left, left),
    inset: (x: 2.5pt, y: 2.35pt),
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else { 0pt },
      bottom: if y == 0 { 0.65pt + ieee_navy } else { 0.32pt + ieee_border },
    ),
    fill: (_, y) => if y == 0 { rgb("#eef2f7") } else if calc.odd(y) { rgb("#fafcff") } else { white },
    table.header(
      [ID], [Algorithm (Year)], [Time Complexity], [Space], [Det.], [Geometry / Bias], [Noise / Outlier], [Auto-$k$], [High-$d$], [Primary Domain],
    ),
    [*176*], [SOM / Kohonen (1982)], [$O(E n M d)$], [$O(M d)$], [Seed], [2D lattice manifold map], [Smooths local outliers], [No ($M$ grid)], [Good (vis.)], [Nonlinear vis., industrial],
    [*177*], [Growing Neural Gas (1995)], [$O(T M_"max" d)$], [$O(M_"max" d + |E|)$], [Order], [Adaptive topology graph], [Age-pruning drops noise], [Yes (comp.)], [Good (manif.)], [Online manifold / robotics],
    [*178*], [ART / Fuzzy ART (1987)], [$O(T k d)$], [$O(k d)$], [Order], [Vigilance hyper-boxes], [Spawns singleton boxes], [Yes ($rho$)], [Fair (comp.)], [Continual non-forgetting],
    [*179*], [Neural Gas (1993)], [$O(T M (d + log M))$], [$O(M d)$], [Stoch.], [Rank-ordered soft VQ], [Moderate (rank decay)], [No ($M$)], [Good (VQ)], [Vector quant., topology],
    [*180*], [Growing Cell Struct. (1994)], [$O(T M_"max" d)$], [$O(M_"max" d)$], [Order], [Fixed-$k$ simplicial mesh], [Prunes low-density cells], [Yes (grow)], [Fair (simp.)], [2D/3D surface modeling],
    [*181*], [DEC (2016)], [$O(E n (W_"enc" + k p))$], [$O(|theta| + B d)$], [Stoch.], [Student-$t$ latent KL], [Latent drift risk (no AE)], [No], [Strong (deep)], [Image / text deep embed.],
    [*182*], [IDEC (2017)], [$O(E n (W_"AE" + k p))$], [$O(|theta| + B d)$], [Stoch.], [Recon.-preserved latent], [AE loss resists collapse], [No], [Strong (deep)], [Tabular & vision deep clus.],
    [*183*], [DCN (2017)], [$O(E n (W_"AE" + k p))$], [$O(|theta| + B d)$], [Stoch.], [Joint AE + $k$-means latent], [AE regularization], [No], [Strong (deep)], [Text, audio, image latent],
    [*184*], [VaDE (2017)], [$O(E n (P + k p))$], [$O(|theta| + B d)$], [Stoch.], [Generative VAE + GMM], [Probabilistic uncertainty], [No], [Strong (gen.)], [Generative image / scRNA],
    [*185*], [DEPICT (2017)], [$O(E n (P + k p))$], [$O(|theta| + B d)$], [Stoch.], [Denoising AE + bal. KL], [Denoising + uniform prior], [No], [Strong (deep)], [Balanced vision datasets],
    [*186*], [JULE (2016)], [$O(T (n P + n^2 d))$], [$O(|theta| + n k_s)$], [Stoch.], [Recurrent CNN + agglom.], [Local $k$NN affinity], [No (stop)], [Strong (CNN)], [Face & object image sets],
    [*187*], [DAC (2017)], [$O(E (n P + n B k))$], [$O(|theta| + B^2)$], [Stoch.], [Pairwise cosine simplex], [Adaptive margin filters], [No], [Strong (CNN)], [Self-paced image clustering],
    [*188*], [IIC (2019)], [$O(E n (C_"net" + k^2))$], [$O(|theta| + B k + k^2)$], [Stoch.], [Augment-invariant MI], [Max-entropy anti-collapse], [No], [Strong (MI)], [Unsupervised vision / seg.],
    [*189*], [SCAN (2020)], [$O(n k log n + E n C_"net")$], [$O(|theta| + n k)$], [Stoch.], [SSL semantic neighbor], [Self-labeling filters bad NN], [No], [Strong (SSL)], [SOTA ImageNet / CIFAR],
    [*190*], [SwAV (2020)], [$O(E n (C_"net" + M k r))$], [$O(|theta| + B_Q r)$], [Stoch.], [Sinkhorn OT prototypes], [Equipartition anti-collapse], [No], [Strong (OT)], [Web-scale SSL vision],
    [*191*], [DeepCluster (2018)], [$O(E (n C_"net" + I n k r'))$], [$O(|theta| + n r')$], [Stoch.], [Alternating CNN + $k$-means], [Empty-cluster split guard], [No], [Strong (CNN)], [Large-scale visual pretrain.],
    [*192*], [Contrastive Clust. / CC (2021)], [$O(E n (C_"enc" + B m + k^2))$], [$O(|theta| + B m)$], [Stoch.], [Dual inst./clus. InfoNCE], [Entropy + contrastive reg.], [No], [Strong (InfoNCE)], [End-to-end vision & text],
    [*193*], [CluStream (2003)], [$O(n q d + I q k d)$], [$O(q d)$ RAM], [Order], [Pyramidal $L_2$ micro-clus.], [Aging deletes old outliers], [No (offline)], [Fair (stream)], [Network intrusion, IoT],
    [*194*], [DenStream (2006)], [$O(n (|P| + |O|) d)$], [$O((|P| + |O|) d)$], [Order], [Fading density micro-clus.], [Explicit $o$-micro-clusters], [Yes ($epsilon$)], [Poor ($epsilon$-ball)], [Non-convex evolving streams],
    [*195*], [StreamKM++ (2012)], [$O(n m d + I m k d)$], [$O(m d log(n \/ m))$], [Seed], [Merge-reduce $L_2$ coreset], [Fragile ($L_2$ coreset)], [No], [Fair (stream)], [Guaranteed stream $k$-means],
    [*196*], [ClusTree (2011)], [$O(n d M log_m K)$], [$O(K d)$ tree], [Order], [Anytime decayed CF-tree], [Exponential time decay], [No (offline)], [Fair (stream)], [High-velocity burst streams],
    [*197*], [D-Stream (2007)], [$O(n d (1 + |cal(G)| \/ t_g))$], [$O(|cal(G)_"act"| d)$], [Order], [Decayed density grid], [Sporadic-grid pruning], [Yes (grid)], [Poor (grid)], [Low-$d$ real-time telemetry],
    [*198*], [EVO-KMeans / GKA (1999)], [$O(G P n k d)$], [$O(P n + P k d + n d)$], [Stoch.], [Global $L_2$ Voronoi (GA)], [Fragile ($L_2$ fitness)], [No], [Poor ($L_2$)], [Small global-opt. bench.],
    [*199*], [PSO-Clustering (2003)], [$O(T S n k d)$], [$O(S k d + n d)$], [Stoch.], [Swarm-searched $L_2$], [Fragile ($L_2$ fitness)], [No], [Poor ($L_2$)], [Continuous global centroid],
    [*200*], [ACO-Clustering (2004)], [$O(T n (R (k + d) + L k d))$], [$O(n (d + k + L))$], [Stoch.], [Pheromone $L_2$ search], [Fragile ($L_2$ fitness)], [No], [Poor ($L_2$)], [Combinatorial global opt.],
  )
})

// =============================================================================
// RESTORE TWO-COLUMN IEEE MODE CLEANLY FOR CONCLUDING SYNTHESIS
// =============================================================================
#set page(columns: 2)

#sec("VI", "Cross-Category Pareto Archetypes and Bridge to Chapters XIII–XIV")

By Corollary 1, no universal "best" clustering algorithm exists across all ten dimensions of Tables IV–XI. Instead, practical mastery of the 200-algorithm compendium reduces to selecting the *Pareto-optimal champion* aligned with the user's binding constraint regime. Synthesizing all eight master matrices yields eight canonical Pareto archetypes spanning the entire theoretical landscape:

1. *Ultra-Fast Streaming & Web-Scale Tabular ($n > 10^7$, fixed $k$):* Mini-Batch K-Means (Alg. 3 @sculley2010webscale), K-Means|| (Alg. 23 @bahmani2012scalable), BICO (Alg. 78 @fichtenberger2013bico), and StreamKM++ (Alg. 195 @ackermann2012streamkmpp) dominate the Tier 1 / Sublinear-RAM frontier; when $n$ fits in linear coordinate RAM, K-Means++ (Alg. 2) paired with Elkan/Hamerly bounds (Algs. 13–14) delivers exact $L_2$ partitioning with Theorem 2's $O(log k)$ guarantee.
2. *Noisy, Non-Convex, Multi-Density Spatial & Exploratory Data ($d <= 15$):* HDBSCAN (Alg. 32 @campello2013hdbscan) dominates fixed-radius DBSCAN (Alg. 31) by eliminating the global $epsilon$ scale trap while retaining $O(n log n)$ indexed speed, strict determinism, automatic $k$ via Excess-of-Mass, and native noise excision $C_(-1)$. For strictly real-time low-$d$ vision/LiDAR, SLIC (Alg. 166) and DBSCAN-3D (Alg. 168) achieve linear/log-linear throughput.
3. *Arbitrary Metric Spaces & Outlier-Contaminated Tabular Data:* FasterPAM (Alg. 24 @schubert2019fasterpam) removes the $O(k)$ slowdown of classical PAM (Alg. 4) for $n <= 10^5$, CLARA-MED (Alg. 25) scales medoids linearly beyond $n > 10^5$, Trimmed K-Means (Alg. 26 @cuesta1997trimmed) enforces impartial $alpha$-trimming, and Student-$t$ Mixtures (Alg. 88 @peel2000robust) provide smooth ellipsoidal generative densities with bounded influence.
4. *Multi-Resolution Hierarchies & Extreme Entity Resolution:* Ward's Method (Alg. 59 @ward1963hierarchical) and Genie (Alg. 69 @gagolewski2016genie) lead deterministic $O(n^2)$ batch dendrograms (compact variance vs. balanced non-convex shapes), Minimax Linkage (Alg. 73 @bien2011minimax) adds interpretable real-point prototypes, and GRINCH / SCC / FINCH (Algs. 77, 79–80) push hierarchical trees to million-scale high-$d$ embeddings.
5. *Probabilistic Model Selection & High-$d$ Text / Count Corpora:* MCLUST (Alg. 93 @fraley2002mclust) and DPGMM (Alg. 82) dominate continuous Gaussian model selection via BIC and Dirichlet-process priors; on sparse high-$d$ bags-of-words and embeddings, Spherical K-Means (Alg. 12), vMF Mixtures (Alg. 89 @banerjee2005vmf), LDA / HDP (Algs. 83, 91), and DMM / BTM (Algs. 97, 101) defeat Euclidean concentration cleanly.
6. *Sparse Complex Networks & Overlapping Communities:* Leiden (Alg. 108 @traag2019leiden) and Infomap (Alg. 113 @rosvall2008infomap) dominate billion-edge disjoint community detection in $O(m)$–$O(m log n)$ time and $O(n + m)$ CSR memory; DC-SBM (Alg. 131 @karrer2011dcsbm) provides rigorous generative inference on scale-free networks; BigCLAM (Alg. 125 @yang2013bigclam) and SLPA (Alg. 126) scale overlapping covers to millions of vertices.
7. *High-Dimensional Subspaces, Two-Mode Genomics & Temporal Series:* PROCLUS / ORCLUS (Algs. 136, 138), P3C (Alg. 148), and Sparse Subspace Clustering (Alg. 150 @elhamifar2013ssc) recover axis-parallel and oblique subspaces buried in ambient noise; FABIA (Alg. 161 @hochreiter2010fabia) leads multiplicative gene biclustering; k-Shape (Alg. 162 @paparrizos2015kshape) and ROCKET-Clustering (Alg. 165 @dempster2020rocket) dominate scalable time-series partitioning.
8. *Raw Perceptual Modalities (Vision, Audio, Multi-Omics):* Self-supervised deep engines—SCAN (Alg. 189 @vangansbeke2020scan), SwAV (Alg. 190 @caron2020swav), Contrastive Clustering (Alg. 192 @li2021contrastive), and generative VaDE (Alg. 184 @jiang2017vade)—jointly learn nonlinear semantic manifolds $Phi_omega$ while preventing representation collapse via Sinkhorn optimal transport, neighbor consistency, or variational priors.

#remark(title: "From Theoretical Reference Matrices to Empirical Benchmarks (Chapter XIII) & Decision Architecture (Chapter XIV)")[
  Tables IV–XI codify the *asymptotic and theoretical* behavior of all 200 clustering algorithms. However, hidden constants, cache locality, distance concentration at moderate $d = 20$–$100$, and real-world cluster overlap routinely alter empirical rankings within each complexity tier. Comparison Chapter XIII therefore puts these theoretical predictions to direct quantitative test through an 800-run reproducible benchmark suite executed on the native Rust `ClusterOverview` engine across spherical, non-convex multi-shape, noisy high-dimensional, and real tabular datasets—before Comparison Chapter XIV synthesizes both theory (Chapter XII) and empirical measurements (Chapter XIII) into a formal ten-vertical engineering decision architecture.
]
