#import "template/ieee.typ": *

// =============================================================================
// IEEE Transactions Monograph Title Banner (Full-Width Top Span)
// =============================================================================
#place(
  top + center,
  scope: "parent",
  float: true,
  block(
    width: 100%,
    above: 0pt,
    below: 14pt,
    {
      set align(center)
      set par(first-line-indent: 0pt, justify: false)
      text(size: 8.5pt, weight: "bold", tracking: 1.2pt, fill: ieee_slate)[
        IEEE TRANSACTIONS MONOGRAPH SERIES  •  UNSUPERVISED LEARNING & MATHEMATICAL CLUSTERING
      ]
      v(6pt)
      line(length: 100%, stroke: 1.4pt + ieee_navy)
      v(10pt)
      text(size: 20pt, weight: "bold", fill: ieee_navy)[
        Mathematical Foundations, Formal Proofs, and Applications of 200 Clustering Algorithms
      ]
      v(5pt)
      text(size: 11pt, style: "italic", fill: ieee_slate)[
        Original Constructive Proofs, Verified Classical Derivations, and Cross-Paradigm Theoretical & Empirical Comparisons
      ]
      v(10pt)
      text(size: 12pt, weight: "bold", fill: rgb("#0f172a"))[
        Created by Gemini, AI made by Google
      ]
      v(2pt)
      text(size: 9pt, fill: ieee_slate)[
        ClusterOverview Mathematical Reference Compendium  •  Edition 1.0 (2026)  •  Companion to the ClusterOverview Native Rust Platform
      ]
      v(8pt)
      line(length: 100%, stroke: 0.7pt + ieee_border)
      v(6pt)
    },
  ),
)

// =============================================================================
// IEEE Abstract & Index Terms
// =============================================================================
#set par(first-line-indent: 0pt)
#text(size: 8.8pt)[
  *Abstract—*This monograph provides a unified, mathematically rigorous IEEE-formatted reference covering all *200 clustering algorithms* spanning eight foundational paradigms: (I) centroid, medoid, and partitioning methods; (II) density-based and mode-seeking methods; (III) hierarchical and tree-based methods; (IV) statistical distribution, mixture, and topic models; (V) graph, spectral, and network community detection; (VI) grid, subspace, and projected clustering; (VII) fuzzy, biclustering, time-series, trajectory, ensemble, and multi-view methods; and (VIII) neural, deep representation, data-stream, and evolutionary metaheuristic clustering. Each algorithm is treated in a dedicated four-page IEEE technical study containing a formal variational or combinatorial specification, pseudocode, explicit theorems and step-by-step proofs—transparently distinguishing *original constructive proofs and unified derivations* from *reproduced or adapted classical theorems with verified scholarly citations*—alongside exact complexity, invariance, breakdown-point, and failure-mode analyses, and concrete engineering and scientific applications. Part IX synthesizes the entire 200-algorithm landscape through fourteen dedicated comparative chapters proving formal cross-family equivalences and asymptotic reductions, axiomatic impossibility results, master 200-algorithm theoretical property matrices, and empirical multi-geometry benchmarks executed via the native Rust `ClusterOverview` engine.

  #v(4pt)
  *Index Terms—*Clustering algorithms, mathematical proofs, Bregman divergence, spectral graph theory, density level sets, ultrametric hierarchies, expectation–maximization, subspace clustering, deep clustering, data streams, comparative analysis, IEEE monograph.
]

#set par(first-line-indent: 1em)

#sec("I", "Authorship, Provenance & Mathematical Protocol")

This monograph was authored and typeset end-to-end as *Created by Gemini, AI made by Google*, designed as the formal mathematical companion to the native Rust `ClusterOverview` platform. Because clustering spans six decades of disparate literatures—information theory, multivariate statistics, algebraic graph theory, computational geometry, signal processing, bioinformatics, and deep representation learning—practitioners and researchers rarely encounter a single work that treats all major methods under one consistent notation and proof standard.

#subsec("A", "Proof Provenance Discipline")

Every formal statement (`Theorem`, `Lemma`, `Proposition`, and `Corollary`) in this compendium carries an explicit provenance tag directly beside its title:

- *Original Proof / Constructive Derivation:* Denotes a proof constructed directly from first principles for this monograph—including unified variational/Bregman reformulations, Lyapunov monotonicity arguments, explicit adversarial counterexamples, exact sensitivity/invariance proofs, finite-sample complexity bounds, and cross-family reduction proofs.
- *Classical / Cited Proof (Ref. [X]):* Denotes an established foundational theorem whose proof is reproduced or rigorously adapted from the cited primary literature (e.g., Lloyd @lloyd1982least, Arthur and Vassilvitskii @arthur2007kmeanspp, Ester et al. @ester1996dbscan, von Luxburg @vonluxburg2007tutorial, Dempster et al. @dempster1977em, Kleinberg @kleinberg2002impossibility). Every bibliography entry corresponds to a verified, published scholarly reference.

#subsec("B", "Standing Mathematical Notation")

Throughout all 200 algorithm entries and 14 comparative chapters, we enforce uniform notation unless locally overridden:

#figure(
  table(
    columns: (1.02fr, 1.98fr),
    align: (left, left),
    [Symbol], [Meaning and Domain],
    [$cal(X) = {x_1, dots, x_n} subset.eq bb(R)^d$], [Finite dataset of $n$ observations in $d$-dimensional Euclidean (or general metric) space $(cal(X), delta)$.],
    [$cal(C) = {C_1, dots, C_k}$], [Hard partition of ${1, dots, n}$ into $k$ disjoint clusters (with optional noise set $C_(-1)$).],
    [$mu_j in bb(R)^d, m_j in cal(X)$], [Centroid (mean prototype) and medoid (exemplar) of cluster $C_j$, with size $n_j = |C_j|$.],
    [$D_phi (x || y)$], [Bregman divergence generated by strictly convex differentiable function $phi$.],
    [$G = (V, E, W)$], [Weighted undirected similarity graph with adjacency $W$, degree matrix $D = "diag"(d_i)$, and Laplacians $L = D - W$, $L_"sym"$, $L_"rw"$.],
    [$U in [0,1]^(n times k)$], [Soft membership / responsibility / posterior assignment matrix with row sums $sum_j u_(i j) = 1$ (relaxed in possibilistic models).],
    [$hat(f)_h (x)$], [Kernel density estimator with bandwidth $h > 0$ and radial kernel $K$.],
    [$J(dot), cal(L)(dot), Q(dot)$], [Clustering objective functional, log-likelihood / ELBO, and Newman–Girvan modularity.],
    [$"ARI", "Sil", "DBI", "CHI"$], [Adjusted Rand Index @hubert1985comparing, Silhouette @rousseeuw1987silhouettes, Davies–Bouldin @davies1979cluster, and Calinski–Harabasz @calinski1974dendrite.],
  ),
  caption: [Standing mathematical notation used uniformly across Parts I–IX.],
)

#sec("II", "Structure of the Compendium")

The monograph is organized into nine parts:
1. *Part I (Algorithms 1–30):* Centroid, Medoid, and Partitioning Algorithms.
2. *Part II (Algorithms 31–55):* Density-Based and Mode-Seeking Algorithms.
3. *Part III (Algorithms 56–80):* Hierarchical and Tree-Based Algorithms.
4. *Part IV (Algorithms 81–102):* Distribution, Statistical Mixture, and Topic Models.
5. *Part V (Algorithms 103–132):* Graph, Spectral, and Network Community Detection.
6. *Part VI (Algorithms 133–150):* Grid, Subspace, and Projected Clustering.
7. *Part VII (Algorithms 151–175):* Fuzzy, Biclustering, Time-Series, Trajectory, Ensemble, and Multi-View Methods.
8. *Part VIII (Algorithms 176–200):* Neural, Deep Embedding, Data-Stream, and Evolutionary Clustering.
9. *Part IX (Comparative Chapters I–XIV):* Grand Unification Theory, Formal Cross-Family Equivalence Proofs, Intra-Category Head-to-Head Comparisons, Axiomatic Impossibility Results, Complete 200-Algorithm Master Property Matrices, Real Multi-Geometry Empirical Benchmarks, and Domain Applications Handbook.

#pagebreak()

// =============================================================================
// Full Monograph Table of Contents
// =============================================================================
#place(
  top + center,
  scope: "parent",
  float: true,
  block(
    width: 100%,
    below: 8pt,
    {
      set align(center)
      text(size: 13pt, weight: "bold", fill: ieee_navy)[TABLE OF CONTENTS]
      v(3pt)
      line(length: 100%, stroke: 0.9pt + ieee_navy)
    },
  ),
)

#outline(
  title: none,
  depth: 2,
  indent: 1em,
)
