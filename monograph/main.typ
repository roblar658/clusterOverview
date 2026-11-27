#import "template/ieee.typ": *
#show: ieee_doc.with()

// Frontmatter: IEEE Title, Abstract, Notation, Provenance & TOC
#include "frontmatter.typ"

#part_opener(
  roman: "I",
  title: "Centroid, Medoid, and Partitioning Algorithms",
  alg_range: "Algorithms 1–30",
  summary: [Part I examines thirty prototype-based partitioning algorithms spanning Euclidean, Manhattan, spherical, Bregman, kernelized, trimmed, weighted, constrained, deterministic, and parallelizable vector-quantization regimes.],
  core_principle: [Every method in Part I minimizes an empirical distortion functional over a discrete Voronoi or exemplar partition, unified by alternating projection onto prototype manifolds and convex or Bregman cell decompositions.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#1.* K-Means (Lloyd's algoritme) \ *\#2.* K-Means++ \ *\#3.* Mini-Batch K-Means \ *\#4.* K-Medoids (PAM) \ *\#5.* CLARA \ *\#6.* CLARANS \ *\#7.* K-Modes \ *\#8.* K-Prototypes \ *\#9.* Bisecting K-Means \ *\#10.* X-Means \ *\#11.* G-Means \ *\#12.* Spherical K-Means \ *\#13.* Elkan's K-Means \ *\#14.* Hamerly's K-Means \ *\#15.* K-Medians], [*\#16.* Kernel K-Means \ *\#17.* Global K-Means \ *\#18.* ISODATA \ *\#19.* Hartigan-Wong K-Means \ *\#20.* COP-KMeans (Constrained K-Means) \ *\#21.* MacQueen's K-Means \ *\#22.* Forgy's K-Means \ *\#23.* K-Means|| (Scalable K-Means++) \ *\#24.* FastPAM / FasterPAM \ *\#25.* CLARA-MED \ *\#26.* Trimmed K-Means \ *\#27.* Power K-Means \ *\#28.* W-K-Means (Weighted K-Means) \ *\#29.* Minkowski Weighted K-Means \ *\#30.* LBG-algoritmen (Linde-Buzo-Gray)]),
)
#include "chapters/alg_001.typ"
#include "chapters/alg_002.typ"
#include "chapters/alg_003.typ"
#include "chapters/alg_004.typ"
#include "chapters/alg_005.typ"
#include "chapters/alg_006.typ"
#include "chapters/alg_007.typ"
#include "chapters/alg_008.typ"
#include "chapters/alg_009.typ"
#include "chapters/alg_010.typ"
#include "chapters/alg_011.typ"
#include "chapters/alg_012.typ"
#include "chapters/alg_013.typ"
#include "chapters/alg_014.typ"
#include "chapters/alg_015.typ"
#include "chapters/alg_016.typ"
#include "chapters/alg_017.typ"
#include "chapters/alg_018.typ"
#include "chapters/alg_019.typ"
#include "chapters/alg_020.typ"
#include "chapters/alg_021.typ"
#include "chapters/alg_022.typ"
#include "chapters/alg_023.typ"
#include "chapters/alg_024.typ"
#include "chapters/alg_025.typ"
#include "chapters/alg_026.typ"
#include "chapters/alg_027.typ"
#include "chapters/alg_028.typ"
#include "chapters/alg_029.typ"
#include "chapters/alg_030.typ"

#part_opener(
  roman: "II",
  title: "Density-Based and Mode-Seeking Algorithms",
  alg_range: "Algorithms 31–55",
  summary: [Part II investigates twenty-five non-parametric density-connectivity, level-set, kernel gradient-ascent, and density-peak algorithms capable of recovering arbitrary non-convex manifolds while isolating background noise.],
  core_principle: [Clusters are formulated either as connected components of upper density level sets {x : f(x) >= lambda} or as basins of attraction of local maxima under continuous or discrete kernel density gradient flows.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#31.* DBSCAN \ *\#32.* HDBSCAN \ *\#33.* OPTICS \ *\#34.* Mean Shift \ *\#35.* DENCLUE \ *\#36.* SNN (Shared Nearest Neighbor) \ *\#37.* PDBSCAN \ *\#38.* DBCLASD \ *\#39.* GDBSCAN \ *\#40.* ST-DBSCAN \ *\#41.* Quickshift \ *\#42.* Medoidshift \ *\#43.* DENCLUE 2.0], [*\#44.* K-DBSCAN \ *\#45.* RNN-DBSCAN \ *\#46.* DeBaCl \ *\#47.* DBSCAN++ \ *\#48.* GriDBSCAN \ *\#49.* MR-DBSCAN \ *\#50.* VDBSCAN \ *\#51.* DPC (Density Peaks Clustering) \ *\#52.* SNN-DPC \ *\#53.* OPTICS-OF \ *\#54.* Blurring Mean Shift \ *\#55.* Mean Shift++]),
)
#include "chapters/alg_031.typ"
#include "chapters/alg_032.typ"
#include "chapters/alg_033.typ"
#include "chapters/alg_034.typ"
#include "chapters/alg_035.typ"
#include "chapters/alg_036.typ"
#include "chapters/alg_037.typ"
#include "chapters/alg_038.typ"
#include "chapters/alg_039.typ"
#include "chapters/alg_040.typ"
#include "chapters/alg_041.typ"
#include "chapters/alg_042.typ"
#include "chapters/alg_043.typ"
#include "chapters/alg_044.typ"
#include "chapters/alg_045.typ"
#include "chapters/alg_046.typ"
#include "chapters/alg_047.typ"
#include "chapters/alg_048.typ"
#include "chapters/alg_049.typ"
#include "chapters/alg_050.typ"
#include "chapters/alg_051.typ"
#include "chapters/alg_052.typ"
#include "chapters/alg_053.typ"
#include "chapters/alg_054.typ"
#include "chapters/alg_055.typ"

#part_opener(
  roman: "III",
  title: "Hierarchical and Tree-Based Algorithms",
  alg_range: "Algorithms 56–80",
  summary: [Part III analyzes twenty-five agglomerative, divisive, metric-tree, phylogenetic, and extreme-scale hierarchical clustering methods constructing multi-resolution nested taxonomies.],
  core_principle: [Hierarchical clustering maps finite dissimilarity spaces into ultrametric trees (dendrograms) governed by the Lance–Williams recursion, minimum spanning tree equivalence, or hierarchical coreset summarization.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#56.* Single Linkage (SLINK) \ *\#57.* Complete Linkage (CLINK) \ *\#58.* Average Linkage (UPGMA) \ *\#59.* Ward's Method \ *\#60.* WPGMA \ *\#61.* Centroid Linkage (UPGMC) \ *\#62.* Median Linkage (WPGMC) \ *\#63.* DIANA \ *\#64.* BIRCH \ *\#65.* CURE \ *\#66.* CHAMELEON \ *\#67.* ROCK \ *\#68.* MONA], [*\#69.* Genie \ *\#70.* HCA-DBSCAN \ *\#71.* BUBBLE \ *\#72.* Lance-Williams algoritme \ *\#73.* Minimax Linkage \ *\#74.* Flexible Beta Linkage \ *\#75.* Neighbor-Joining (NJ) \ *\#76.* PERCH \ *\#77.* GRINCH \ *\#78.* BICO \ *\#79.* SCC (Sub-Cluster Component) \ *\#80.* Finch]),
)
#include "chapters/alg_056.typ"
#include "chapters/alg_057.typ"
#include "chapters/alg_058.typ"
#include "chapters/alg_059.typ"
#include "chapters/alg_060.typ"
#include "chapters/alg_061.typ"
#include "chapters/alg_062.typ"
#include "chapters/alg_063.typ"
#include "chapters/alg_064.typ"
#include "chapters/alg_065.typ"
#include "chapters/alg_066.typ"
#include "chapters/alg_067.typ"
#include "chapters/alg_068.typ"
#include "chapters/alg_069.typ"
#include "chapters/alg_070.typ"
#include "chapters/alg_071.typ"
#include "chapters/alg_072.typ"
#include "chapters/alg_073.typ"
#include "chapters/alg_074.typ"
#include "chapters/alg_075.typ"
#include "chapters/alg_076.typ"
#include "chapters/alg_077.typ"
#include "chapters/alg_078.typ"
#include "chapters/alg_079.typ"
#include "chapters/alg_080.typ"

#part_opener(
  roman: "IV",
  title: "Distribution, Statistical Mixture, and Topic Models",
  alg_range: "Algorithms 81–102",
  summary: [Part IV develops twenty-two generative probabilistic mixture models, Bayesian nonparametric processes, conceptual utility trees, and latent topic/factorization models.],
  core_principle: [Observations arise from latent-variable mixture densities maximized via Expectation–Maximization, variational lower bounds (ELBO), collapsed Gibbs sampling, Bayesian Information Criterion, or Minimum Message Length.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#81.* Gaussian Mixture Models (GMM) \ *\#82.* Bayesian GMM (DPGMM) \ *\#83.* Latent Dirichlet Allocation (LDA) \ *\#84.* pLSA \ *\#85.* COBWEB \ *\#86.* CLASSIT \ *\#87.* AutoClass \ *\#88.* Student's t-Mixture Models (TMM) \ *\#89.* von Mises-Fisher Mixture Models (vMF) \ *\#90.* Latent Class Analysis (LCA) \ *\#91.* Hierarchical Dirichlet Process (HDP)], [*\#92.* NMF (Non-Negative Matrix Factorization) \ *\#93.* MCLUST \ *\#94.* Mixmod \ *\#95.* Poisson Mixture Models \ *\#96.* Bernoulli Mixture Models \ *\#97.* Dirichlet Multinomial Mixture (DMM) \ *\#98.* Skew-Normal Mixture Models \ *\#99.* Infinite Gaussian Mixture Model (IGMM) \ *\#100.* Correlated Topic Model (CTM) \ *\#101.* Biterm Topic Model (BTM) \ *\#102.* SNOB]),
)
#include "chapters/alg_081.typ"
#include "chapters/alg_082.typ"
#include "chapters/alg_083.typ"
#include "chapters/alg_084.typ"
#include "chapters/alg_085.typ"
#include "chapters/alg_086.typ"
#include "chapters/alg_087.typ"
#include "chapters/alg_088.typ"
#include "chapters/alg_089.typ"
#include "chapters/alg_090.typ"
#include "chapters/alg_091.typ"
#include "chapters/alg_092.typ"
#include "chapters/alg_093.typ"
#include "chapters/alg_094.typ"
#include "chapters/alg_095.typ"
#include "chapters/alg_096.typ"
#include "chapters/alg_097.typ"
#include "chapters/alg_098.typ"
#include "chapters/alg_099.typ"
#include "chapters/alg_100.typ"
#include "chapters/alg_101.typ"
#include "chapters/alg_102.typ"

#part_opener(
  roman: "V",
  title: "Graph, Spectral, and Network Community Detection",
  alg_range: "Algorithms 103–132",
  summary: [Part V covers thirty spectral graph cut, modularity optimization, random-walk flow, message-passing, overlapping community, multilevel partitioning, and stochastic blockmodel algorithms.],
  core_principle: [Partitions optimize discrete graph cut ratios, Newman–Girvan modularity, Markov random-walk description length, or generative block-edge likelihoods via Laplacian eigen-relaxations and local greedy/message updates.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#103.* Spectral Clustering (Ng-Jordan-Weiss) \ *\#104.* Normalized Cuts (Ncut) \ *\#105.* Ratio Cut \ *\#106.* Affinity Propagation \ *\#107.* Louvain Method \ *\#108.* Leiden Algorithm \ *\#109.* Markov Clustering (MCL) \ *\#110.* Label Propagation (LPA) \ *\#111.* Chinese Whispers \ *\#112.* Girvan-Newman \ *\#113.* Infomap \ *\#114.* Walktrap \ *\#115.* Fast Greedy Modularity (CNM) \ *\#116.* Leading Eigenvector (Newman) \ *\#117.* Spinglass Community Detection], [*\#118.* HCS (Highly Connected Subgraphs) \ *\#119.* FLUIDC (Fluid Communities) \ *\#120.* Kernighan-Lin algoritme \ *\#121.* METIS \ *\#122.* Chaco \ *\#123.* Spectral Biclustering (Kluger) \ *\#124.* OSLOM \ *\#125.* BigCLAM \ *\#126.* SLPA (Speaker-Listener Label Propagation) \ *\#127.* CPM (Clique Percolation Method) \ *\#128.* DEMON \ *\#129.* Fast Unfolding (Blondel) \ *\#130.* SBM (Stochastic Block Model) \ *\#131.* Degree-Corrected SBM \ *\#132.* Power Iteration Clustering (PIC)]),
)
#include "chapters/alg_103.typ"
#include "chapters/alg_104.typ"
#include "chapters/alg_105.typ"
#include "chapters/alg_106.typ"
#include "chapters/alg_107.typ"
#include "chapters/alg_108.typ"
#include "chapters/alg_109.typ"
#include "chapters/alg_110.typ"
#include "chapters/alg_111.typ"
#include "chapters/alg_112.typ"
#include "chapters/alg_113.typ"
#include "chapters/alg_114.typ"
#include "chapters/alg_115.typ"
#include "chapters/alg_116.typ"
#include "chapters/alg_117.typ"
#include "chapters/alg_118.typ"
#include "chapters/alg_119.typ"
#include "chapters/alg_120.typ"
#include "chapters/alg_121.typ"
#include "chapters/alg_122.typ"
#include "chapters/alg_123.typ"
#include "chapters/alg_124.typ"
#include "chapters/alg_125.typ"
#include "chapters/alg_126.typ"
#include "chapters/alg_127.typ"
#include "chapters/alg_128.typ"
#include "chapters/alg_129.typ"
#include "chapters/alg_130.typ"
#include "chapters/alg_131.typ"
#include "chapters/alg_132.typ"

#part_opener(
  roman: "VI",
  title: "Grid, Subspace, and Projected Clustering Algorithms",
  alg_range: "Algorithms 133–150",
  summary: [Part VI treats eighteen spatial grid quantization, wavelet multiresolution, axis-parallel subspace, arbitrarily oriented projected, and sparse self-expressive subspace clustering algorithms designed to defeat the curse of dimensionality.],
  core_principle: [High-dimensional structure is recovered by restricting distortion or density evaluation to cluster-specific low-dimensional subspaces via downward-closed pruning, entropy criteria, local PCA, or L1 sparse representation.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#133.* STING \ *\#134.* CLIQUE \ *\#135.* WaveCluster \ *\#136.* PROCLUS \ *\#137.* SUBCLU \ *\#138.* ORCLUS \ *\#139.* DOC \ *\#140.* MAFIA \ *\#141.* OptiGrid], [*\#142.* BANG-Clustering \ *\#143.* ENCLUS \ *\#144.* FINDIT \ *\#145.* COSA (Clustering Objects on Subsets of Attributes) \ *\#146.* FIRES \ *\#147.* INSCY \ *\#148.* P3C \ *\#149.* STATPC \ *\#150.* Sparse Subspace Clustering (SSC)]),
)
#include "chapters/alg_133.typ"
#include "chapters/alg_134.typ"
#include "chapters/alg_135.typ"
#include "chapters/alg_136.typ"
#include "chapters/alg_137.typ"
#include "chapters/alg_138.typ"
#include "chapters/alg_139.typ"
#include "chapters/alg_140.typ"
#include "chapters/alg_141.typ"
#include "chapters/alg_142.typ"
#include "chapters/alg_143.typ"
#include "chapters/alg_144.typ"
#include "chapters/alg_145.typ"
#include "chapters/alg_146.typ"
#include "chapters/alg_147.typ"
#include "chapters/alg_148.typ"
#include "chapters/alg_149.typ"
#include "chapters/alg_150.typ"

#part_opener(
  roman: "VII",
  title: "Fuzzy, Biclustering, Time-Series, and Domain-Specific Methods",
  alg_range: "Algorithms 151–175",
  summary: [Part VII addresses twenty-five soft-membership fuzzy/possibilistic partitioners, two-mode biclustering models, elastic/shape-invariant time-series clusterers, superpixel/3D/trajectory algorithms, consensus ensembles, and multi-view methods.],
  core_principle: [Domain invariances—gradual membership, simultaneous row-column coherence, phase/warp invariance, spatial continuity, and cross-view consensus—are encoded directly into generalized dissimilarity functionals and hypergraph objectives.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#151.* Fuzzy C-Means (FCM) \ *\#152.* Possibilistic C-Means (PCM) \ *\#153.* Gustafson-Kessel (GK) \ *\#154.* Gath-Geva (GG) \ *\#155.* Kernel Fuzzy C-Means (KFCM) \ *\#156.* Spatial Fuzzy C-Means (sFCM) \ *\#157.* Cheng-Church Biclustering \ *\#158.* PLAID \ *\#159.* Bimax \ *\#160.* ISA (Iterative Signature Algorithm) \ *\#161.* FABIA \ *\#162.* k-Shape \ *\#163.* DBA (DTW Barycenter Averaging)], [*\#164.* Soft-DTW Clustering \ *\#165.* ROCKET-Clustering \ *\#166.* SLIC (Simple Linear Iterative Clustering) \ *\#167.* Watershed Clustering \ *\#168.* DBSCAN-3D / Euclidean Cluster Extraction \ *\#169.* TRACLUS \ *\#170.* Consensus Clustering (Monti) \ *\#171.* CSPA (Cluster-based Similarity Partitioning) \ *\#172.* HGPA (HyperGraph Partitioning Algorithm) \ *\#173.* MCLA (Meta-CLustering Algorithm) \ *\#174.* Multi-View K-Means \ *\#175.* Co-regularized Multi-View Spectral Clustering]),
)
#include "chapters/alg_151.typ"
#include "chapters/alg_152.typ"
#include "chapters/alg_153.typ"
#include "chapters/alg_154.typ"
#include "chapters/alg_155.typ"
#include "chapters/alg_156.typ"
#include "chapters/alg_157.typ"
#include "chapters/alg_158.typ"
#include "chapters/alg_159.typ"
#include "chapters/alg_160.typ"
#include "chapters/alg_161.typ"
#include "chapters/alg_162.typ"
#include "chapters/alg_163.typ"
#include "chapters/alg_164.typ"
#include "chapters/alg_165.typ"
#include "chapters/alg_166.typ"
#include "chapters/alg_167.typ"
#include "chapters/alg_168.typ"
#include "chapters/alg_169.typ"
#include "chapters/alg_170.typ"
#include "chapters/alg_171.typ"
#include "chapters/alg_172.typ"
#include "chapters/alg_173.typ"
#include "chapters/alg_174.typ"
#include "chapters/alg_175.typ"

#part_opener(
  roman: "VIII",
  title: "Neural, Deep, Streaming, and Evolutionary Algorithms",
  alg_range: "Algorithms 176–200",
  summary: [Part VIII studies twenty-five competitive neural maps, deep representation-clustering architectures, single-pass evolving data-stream engines, and population-based stochastic global metaheuristics.],
  core_principle: [Clustering is coupled either with nonlinear representation learning and self-supervised mutual-information/optimal-transport losses, bounded-memory fading micro-cluster coresets, or global stochastic search dynamics.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [*\#176.* SOM (Self-Organizing Maps) \ *\#177.* Growing Neural Gas (GNG) \ *\#178.* ART (Adaptive Resonance Theory) \ *\#179.* Neural Gas \ *\#180.* Growing Cell Structures (GCS) \ *\#181.* DEC (Deep Embedded Clustering) \ *\#182.* IDEC (Improved DEC) \ *\#183.* DCN (Deep Clustering Network) \ *\#184.* VaDE (Variational Deep Embedding) \ *\#185.* DEPICT \ *\#186.* JULE (Joint Unsupervised Learning) \ *\#187.* DAC (Deep Adaptive Clustering) \ *\#188.* IIC (Invariant Information Clustering)], [*\#189.* SCAN (Semantic Clustering by Adopting Nearest neighbors) \ *\#190.* SwAV (Swapping Assignments between Views) \ *\#191.* DeepCluster \ *\#192.* Contrastive Clustering (CC) \ *\#193.* CluStream \ *\#194.* DenStream \ *\#195.* StreamKM++ \ *\#196.* BIRCH-Stream / ClusTree \ *\#197.* D-Stream \ *\#198.* EVO-KMeans / Genetic K-Means (GKA) \ *\#199.* PSO-Clustering (Particle Swarm Optimization) \ *\#200.* ACO-Clustering (Ant Colony Optimization)]),
)
#include "chapters/alg_176.typ"
#include "chapters/alg_177.typ"
#include "chapters/alg_178.typ"
#include "chapters/alg_179.typ"
#include "chapters/alg_180.typ"
#include "chapters/alg_181.typ"
#include "chapters/alg_182.typ"
#include "chapters/alg_183.typ"
#include "chapters/alg_184.typ"
#include "chapters/alg_185.typ"
#include "chapters/alg_186.typ"
#include "chapters/alg_187.typ"
#include "chapters/alg_188.typ"
#include "chapters/alg_189.typ"
#include "chapters/alg_190.typ"
#include "chapters/alg_191.typ"
#include "chapters/alg_192.typ"
#include "chapters/alg_193.typ"
#include "chapters/alg_194.typ"
#include "chapters/alg_195.typ"
#include "chapters/alg_196.typ"
#include "chapters/alg_197.typ"
#include "chapters/alg_198.typ"
#include "chapters/alg_199.typ"
#include "chapters/alg_200.typ"

#part_opener(
  roman: "IX",
  title: "Cross-Algorithm Mathematical Unification, Equivalences, and Empirical Comparisons",
  alg_range: "Comparative Chapters I–XIV",
  summary: [Part IX synthesizes the entire 200-algorithm compendium through fourteen dedicated comparative chapters proving cross-family mathematical equivalences, axiomatic impossibility boundaries, intra-category head-to-head analyses, master theoretical reference matrices, real multi-geometry benchmarks, and domain engineering decision architectures.],
  core_principle: [Apparent algorithmic diversity collapses into a tightly connected web of variational limits, spectral-kernel dualities, ultrametric level-set functors, and information-theoretic trade-offs governed by Kleinberg impossibility and high-dimensional measure concentration.],
  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [I. Grand Unified Variational, Geometric & Information-Theoretic Framework \ II. Formal Mathematical Equivalences, Asymptotic Limits & Reductions \ III. Intra-Category Comparison: Centroid, Medoid & Partitioning (1–30) \ IV. Intra-Category Comparison: Density & Mode-Seeking (31–55) \ V. Intra-Category Comparison: Hierarchical & Tree-Based (56–80) \ VI. Intra-Category Comparison: Distribution, Statistical & Topic Models (81–102) \ VII. Intra-Category Comparison: Graph, Spectral & Network (103–132)], [VIII. Intra-Category Comparison: Grid, Subspace & Projection (133–150) \ IX. Intra-Category Comparison: Fuzzy, Biclustering, Time-Series & Domain (151–175) \ X. Intra-Category Comparison: Neural, Deep, Streaming & Evolutionary (176–200) \ XI. Axiomatic Foundations, Impossibility Theorems & Breakdown Robustness \ XII. Complete 200-Algorithm Master Comparative Reference Matrices \ XIII. Empirical Multi-Geometry Benchmark Comparison (ClusterOverview Engine) \ XIV. Practitioner's Engineering Applications Handbook & Selection Architecture]),
)
#include "chapters/cmp_01.typ"
#include "chapters/cmp_02.typ"
#include "chapters/cmp_03.typ"
#include "chapters/cmp_04.typ"
#include "chapters/cmp_05.typ"
#include "chapters/cmp_06.typ"
#include "chapters/cmp_07.typ"
#include "chapters/cmp_08.typ"
#include "chapters/cmp_09.typ"
#include "chapters/cmp_10.typ"
#include "chapters/cmp_11.typ"
#include "chapters/cmp_12.typ"
#include "chapters/cmp_13.typ"
#include "chapters/cmp_14.typ"

#pagebreak(weak: true)
#place(top + center, scope: "parent", float: true, block(width: 100%, fill: ieee_navy, inset: (x: 16pt, y: 12pt), radius: 3pt, below: 10pt, {
  set align(center)
  heading(level: 1, outlined: true, numbering: none)[Complete IEEE Bibliography]
}))
#bibliography("references.bib", title: none, style: "ieee")