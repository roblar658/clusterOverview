import asyncio
import json
import subprocess
from pathlib import Path

MONO = Path(r"C:\dev\clusterOverview\monograph")
CH_DIR = MONO / "chapters"
CH_DIR.mkdir(parents=True, exist_ok=True)

catalog = json.loads((MONO / "data" / "catalog_200.json").read_text(encoding="utf-8"))
bib_keys = json.loads((MONO / "data" / "bib_keys.json").read_text(encoding="utf-8"))
bench = json.loads((MONO / "data" / "benchmark_summary.json").read_text(encoding="utf-8"))
BIB_KEYS_STR = ", ".join(bib_keys)

wave_size = int(args.get("wave_size", 5))


def is_cmp_verified(cid: int, min_p: int, max_p: int):
    f = CH_DIR / f"cmp_{cid:02d}.typ"
    if not f.exists() or f.stat().st_size < 8000:
        return False, 0
    res = subprocess.run(
        ["python", str(MONO / "scripts" / "verify_chapter.py"), f"chapters/cmp_{cid:02d}.typ", str(min_p), str(max_p)],
        cwd=str(MONO),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if res.returncode != 0:
        return False, 0
    try:
        data = json.loads(res.stdout)
        return bool(data.get("ok") and data.get("length_ok") and not data.get("warnings")), int(data.get("content_pages", 0))
    except Exception:
        return False, 0


def cat_list_str(cat_num: int) -> str:
    items = [it for it in catalog if it["cat_num"] == cat_num]
    return "\n".join(f"  - Alg. {it['id']}: {it['name']} ({it['desc_no']})" for it in items)


CMP_SPECS = [
    {
        "cid": 1, "roman": "I", "min_p": 10, "max_p": 18,
        "title": "Grand Unified Variational, Geometric, and Information-Theoretic Framework",
        "scope": "Unifying All Eight Clustering Families (Algorithms 1-200)",
        "instructions": (
            "Write a deep 12-15 page IEEE comparative synthesis unifying ALL 8 clustering families under six master mathematical lenses: "
            "(1) Bregman Divergences & Dual Exponential Families (@banerjee2005clustering), "
            "(2) Variational Free Energy / Helmholtz Functional connecting hard partitioning, EM, fuzzy membership, and deep VAEs, "
            "(3) Topological Level Sets, Persistent Homology & Ultrametric Dendrograms (@carlsson2010characterization), "
            "(4) Spectral Dirichlet Energies, Cheeger Inequalities & Random-Walk Diffusion (@vonluxburg2007tutorial, @chung1997spectral), "
            "(5) Minimum Description Length (MDL), Rate-Distortion & Information Bottleneck (Infomap, SNOB, X-Means BIC, IIC), and "
            "(6) Entropic Optimal Transport & Balanced Assignment (Sinkhorn-Knopp, SwAV). "
            "Include 6-8 formal unifying Theorems/Propositions with full proofs (explicitly marking Original Unified Proofs vs. Cited results), "
            "embed `#wide_figure(caption: [Mathematical unification map across all eight clustering families.])[#image(\"../figures/fig_taxonomy_map.svg\", width: 96%)]`, "
            "and include 3 comprehensive IEEE comparison tables mapping all 8 families across objective geometry, regularizers, and statistical assumptions."
        ),
    },
    {
        "cid": 2, "roman": "II", "min_p": 14, "max_p": 24,
        "title": "Formal Mathematical Equivalences, Asymptotic Limits, and Exact Cross-Family Reductions",
        "scope": "15+ Complete Cross-Family Equivalence & Reduction Proofs Across All 200 Algorithms",
        "instructions": (
            "Write a landmark 16-20 page IEEE chapter containing 14-16 rigorous, complete formal Theorems & Proofs establishing exact mathematical equivalences, asymptotic limits, and reductions connecting disparate families in our 200-algorithm catalog:\n"
            "1. Weighted Kernel K-Means (Alg. 16) <=> Normalized Cut (Alg. 104) & Ratio Cut (Alg. 105) (@dhillon2004kernel) [Complete matrix trace derivation]\n"
            "2. Isotropic GMM EM (Alg. 81) -> Lloyd's K-Means (Alg. 1) as sigma -> 0+ [Original rigorous Laplace-limit proof]\n"
            "3. Fuzzy C-Means (Alg. 151) -> K-Means (Alg. 1) as m -> 1+ and Entropy-Regularized K-Means <=> Softmax EM [Original proof]\n"
            "4. Single Linkage (Alg. 56) <=> Kruskal MST Thresholding <=> Core-Point DBSCAN (Alg. 31) Level-Set Components [Original graph-theoretic proof]\n"
            "5. HDBSCAN (Alg. 32) <=> Single Linkage on Mutual Reachability Metric [Constructive metric-space proof]\n"
            "6. Mean Shift (Alg. 34) <=> EM Monotone Ascent for Gaussian KDE <=> DENCLUE 2.0 (Alg. 43) [Exact fixed-point proof]\n"
            "7. Ward's Method (Alg. 59) <=> Greedy Agglomerative Minimization of K-Means WCSS (Alg. 1) [Exact Lance-Williams variance delta proof]\n"
            "8. KL-NMF (Alg. 92) <=> pLSA (Alg. 84) Maximum Likelihood (@ding2008nmfplsa) [Complete equivalence proof]\n"
            "9. Newman-Girvan Modularity (Louvain/Leiden, Algs. 107-108) <=> Potts Spinglass Hamiltonian (Alg. 117) <=> Uniform-Planting DC-SBM MLE (Alg. 131) [Unified proof]\n"
            "10. Spherical K-Means (Alg. 12) <=> Zero-Temperature Limit (kappa -> inf) of von Mises-Fisher Mixture (Alg. 89) [Original proof]\n"
            "11. Spectral NJW (Alg. 103) <=> Random Walk Commute-Time Embedding (Walktrap Alg. 114, PIC Alg. 132) [Spectral proof]\n"
            "12. SwAV Sinkhorn-Knopp (Alg. 190) <=> Entropic Optimal Transport Balanced K-Means (@cuturi2013sinkhorn, @caron2020swav) [Dual proof]\n"
            "13. BICO (Alg. 78) & StreamKM++ (Alg. 195) Coreset Distortion Preservation for K-Means [Approximation proof]\n"
            "14. Sparse Subspace Clustering (Alg. 150) -> Subspace-Preserving Spectral Cut (Alg. 103) (@elhamifar2013ssc) [Self-expressive proof]\n"
            "15. Deep Clustering (DEC/IDEC/DCN/VaDE, Algs. 181-184) <=> Amortized Variational EM / Alternating Manifold K-Means [Unified proof].\n"
            "Conclude with a full Master Reduction Summary Table mapping every equivalence, exact limit condition, and practical algorithmic implication."
        ),
    },
    {
        "cid": 3, "roman": "III", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Centroid, Medoid, and Partitioning (Algs. 1-30)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 1-30",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 30 algorithms in Category I:\n"
            + cat_list_str(1)
            + "\nOrganize into: (I) Seeding & Approximation Guarantees (Lloyd vs. Forgy vs. MacQueen vs. K-Means++ vs. K-Means|| vs. Global K-Means vs. Power K-Means — with original comparative theorems on approximation ratios & variance), "
            "(II) Exact Acceleration & Geometric Pruning (Elkan vs. Hamerly vs. Hartigan-Wong vs. Mini-Batch — proving exact equivalence of Elkan/Hamerly to Lloyd vs. strict descent advantage of Hartigan-Wong), "
            "(III) Metric Generalizations & Robustness (K-Means vs. K-Medians vs. PAM/FasterPAM vs. CLARA/CLARANS/CLARA-MED vs. Trimmed K-Means vs. Spherical vs. Kernel vs. K-Modes/K-Prototypes — proving breakdown-point hierarchy), "
            "(IV) Automatic Model Selection & Feature Weighting (X-Means vs. G-Means vs. ISODATA vs. Bisecting vs. W-K-Means vs. Minkowski Weighted vs. COP-KMeans vs. LBG), "
            "and (V) Complete 30-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 4, "roman": "IV", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Density and Mode-Seeking (Algs. 31-55)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 31-55",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 25 algorithms in Category II:\n"
            + cat_list_str(2)
            + "\nOrganize into: (I) Density-Reachability & Variable-Density Families (DBSCAN vs. OPTICS vs. HDBSCAN vs. VDBSCAN vs. RNN-DBSCAN vs. SNN vs. GDBSCAN vs. ST-DBSCAN vs. DBCLASD — proving why fixed-eps DBSCAN fails under multi-scale densities while HDBSCAN/OPTICS/RNN-DBSCAN recover nested level sets), "
            "(II) Kernel Mode-Seeking & Gradient Ascent Families (Mean Shift vs. Blurring Mean Shift vs. Mean Shift++ vs. DENCLUE 1.0/2.0 vs. Quickshift vs. Medoidshift — proving step-size, fixed-point, and discrete-medoid tradeoffs), "
            "(III) Density Peaks & Level-Set Trees (DPC vs. SNN-DPC vs. DeBaCl vs. OPTICS-OF — analyzing decision-graph geometry and outlier scoring), "
            "(IV) Parallel, Grid, and Subsampled Scalability (PDBSCAN vs. MR-DBSCAN vs. GriDBSCAN vs. K-DBSCAN vs. DBSCAN++ — proving communication vs. accuracy bounds), "
            "and (V) Complete 25-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 5, "roman": "V", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Hierarchical and Tree-Based (Algs. 56-80)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 56-80",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 25 algorithms in Category III:\n"
            + cat_list_str(3)
            + "\nOrganize into: (I) The Lance-Williams Unification & Monotonicity (SLINK, CLINK, UPGMA, Ward, WPGMA, UPGMC, WPGMC, Flexible Beta — proving exact Milligan/Lance-Williams monotonicity conditions and why Centroid/Median exhibit inversions while Single/Complete/Average/Ward remain strictly ultrametric), "
            "(II) Chaining, Diameter, and Dasgupta Cost Guarantees (Single vs. Complete vs. Average vs. Minimax vs. Genie — with proofs of Gini-controlled anti-chaining and Dasgupta cost bounds @dasgupta2016cost, @moseley2017approximation), "
            "(III) Scalable Multi-Stage & Metric-Space Trees (BIRCH vs. BICO vs. BUBBLE vs. CURE vs. CHAMELEON vs. ROCK vs. HCA-DBSCAN), "
            "(IV) Divisive, Phylogenetic & Extreme-Scale Incremental Hierarchies (DIANA vs. MONA vs. Neighbor-Joining vs. PERCH vs. GRINCH vs. SCC vs. Finch), "
            "and (V) Complete 25-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 6, "roman": "VI", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Distribution, Statistical, and Topic Models (Algs. 81-102)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 81-102",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 22 algorithms in Category IV:\n"
            + cat_list_str(4)
            + "\nOrganize into: (I) Finite Parametric Mixtures & Covariance Parsimony (GMM vs. MCLUST EII-VVV vs. Mixmod vs. Student-t TMM vs. Skew-Normal vs. vMF — proving identifiability, singularity pathologies, and heavy-tail influence-function boundedness of TMM), "
            "(II) Bayesian Nonparametric Mixtures & Model Selection (Bayesian DPGMM vs. IGMM/CRP vs. AutoClass vs. SNOB/MML — contrasting stick-breaking variational Bayes, collapsed Gibbs CRP, BIC, and MML), "
            "(III) Discrete, Count, Binary & Conceptual Mixtures (Poisson Mixture vs. Bernoulli Mixture vs. LCA vs. COBWEB vs. CLASSIT), "
            "(IV) Topic Models & Matrix Factorization (pLSA vs. LDA vs. HDP vs. CTM vs. DMM vs. BTM vs. NMF — proving short-text sparsity collapse of LDA vs. biterm/single-topic remedies of BTM/DMM), "
            "and (V) Complete 22-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 7, "roman": "VII", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Graph, Spectral, and Network Algorithms (Algs. 103-132)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 103-132",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 30 algorithms in Category V:\n"
            + cat_list_str(5)
            + "\nOrganize into: (I) Spectral Cuts & Laplacian Families (NJW vs. Ncut vs. Ratio Cut vs. Spectral Biclustering vs. PIC — proving Cheeger inequality bounds and unnormalized vs. normalized Laplacian consistency @vonluxburg2007tutorial), "
            "(II) Modularity, Resolution Limits & Connectedness (CNM vs. Leading Eigenvector vs. Louvain/Fast Unfolding vs. Leiden vs. Spinglass — proving the Fortunato-Barthélemy resolution limit @fortunato2007resolution and Traag et al.'s disconnected-community defect in Louvain resolved by Leiden @traag2019leiden), "
            "(III) Random-Walk, Flow, Label & Message-Passing Dynamics (MCL vs. Infomap vs. Walktrap vs. LPA vs. Chinese Whispers vs. FLUIDC vs. Affinity Propagation), "
            "(IV) Divisive Cuts, Multilevel Partitioning, Overlapping & Generative Block Models (Girvan-Newman vs. HCS vs. Kernighan-Lin vs. METIS vs. Chaco vs. OSLOM vs. BigCLAM vs. SLPA vs. CPM vs. DEMON vs. SBM vs. DC-SBM — contrasting exact recovery thresholds @abbe2018sbmreview), "
            "and (V) Complete 30-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 8, "roman": "VIII", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Grid, Subspace, and Projection (Algs. 133-150)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 133-150",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 18 algorithms in Category VI:\n"
            + cat_list_str(6)
            + "\nOrganize into: (I) Pure Grid & Multi-Resolution Spatial Partitioning (STING vs. WaveCluster vs. OptiGrid vs. BANG — proving wavelet noise-suppression bounds and optimal low-density cutting-plane guarantees), "
            "(II) Axis-Parallel Bottom-Up Subspace Clustering & Monotonicity (CLIQUE vs. MAFIA vs. SUBCLU vs. ENCLUS vs. FIRES vs. INSCY — proving Apriori downward-closure of density & entropy and redundancy-pruning complexity), "
            "(III) Projected Clustering & Soft Attribute Weighting (PROCLUS vs. DOC vs. FINDIT vs. COSA vs. P3C vs. STATPC — contrasting k-medoid subspace projection, Monte Carlo hypercube sampling, and Poisson statistical significance), "
            "(IV) Arbitrarily Oriented & Sparse Algebraic Subspaces (ORCLUS vs. Sparse Subspace Clustering SSC — proving L1 subspace-detection property vs. local PCA orientation @elhamifar2013ssc, @soltanolkotabi2012geometric), "
            "and (V) Complete 18-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 9, "roman": "IX", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Fuzzy, Biclustering, Time-Series, and Domain Methods (Algs. 151-175)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 151-175",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 25 algorithms in Category VII:\n"
            + cat_list_str(7)
            + "\nOrganize into: (I) Fuzzy & Possibilistic Prototype Families (FCM vs. PCM vs. Gustafson-Kessel vs. Gath-Geva vs. KFCM vs. sFCM — proving why PCM eliminates the sum-to-one noise trap of FCM yet risks coincident clusters, and comparing spherical/Mahalanobis/kernel/spatial metrics), "
            "(II) Biclustering Models for Two-Mode Matrices (Cheng-Church vs. PLAID vs. Bimax vs. ISA vs. FABIA — contrasting constant/coherent MSR, additive layers, binary submatrices, singular-vector signatures, and sparse factor models), "
            "(III) Time-Series Elastic & Shape Invariance (k-Shape vs. DBA vs. Soft-DTW vs. ROCKET-Clustering — proving phase/shift/amplitude vs. nonlinear warp invariances and Fréchet mean smoothness), "
            "(IV) Spatial Vision, 3D Point Clouds, Trajectories, Ensembles & Multi-View (SLIC vs. Watershed vs. DBSCAN-3D vs. TRACLUS vs. Consensus Monti vs. CSPA/HGPA/MCLA vs. Multi-View K-Means vs. Co-regularized Multi-View Spectral — proving consensus information-theoretic optimality @strehl2002cluster), "
            "and (V) Complete 25-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 10, "roman": "X", "min_p": 8, "max_p": 14,
        "title": "Intra-Category Comparative Analysis: Neural, Deep, Streaming, and Evolutionary (Algs. 176-200)",
        "scope": "Head-to-Head Theoretical & Practical Comparison of Algorithms 176-200",
        "instructions": (
            "Write a deep 10-12 page IEEE comparative study contrasting ALL 25 algorithms in Category VIII:\n"
            + cat_list_str(8)
            + "\nOrganize into: (I) Classical Competitive & Topological Neural Networks (SOM vs. Neural Gas vs. GNG vs. GCS vs. ART/Fuzzy ART — proving topological ordering vs. stability-plasticity vigilance guarantees), "
            "(II) Deep Representation & Self-Supervised Clustering (DEC vs. IDEC vs. DCN vs. VaDE vs. DEPICT vs. JULE vs. DAC vs. IIC vs. SCAN vs. SwAV vs. DeepCluster vs. Contrastive Clustering — proving feature-collapse modes of DEC/DeepCluster and how reconstruction loss in IDEC/DCN, ELBO in VaDE, Sinkhorn equipartition in SwAV, and InfoNCE mutual information in IIC/CC prevent trivial collapse), "
            "(III) Single-Pass & Anytime Data-Stream Engines (CluStream vs. DenStream vs. StreamKM++ vs. ClusTree vs. D-Stream — proving exponential decay weight bounds, micro-cluster additivity, and coreset distortion guarantees), "
            "(IV) Global Stochastic Metaheuristics (EVO-KMeans/GKA vs. PSO-Clustering vs. ACO-Clustering — proving Markov-chain global convergence in GKA vs. finite-time premature convergence risks), "
            "and (V) Complete 25-Algorithm Head-to-Head Comparison Table & Decision Guide."
        ),
    },
    {
        "cid": 11, "roman": "XI", "min_p": 10, "max_p": 16,
        "title": "Axiomatic Foundations, Impossibility Theorems, Concentration of Measure, and Robustness",
        "scope": "Cross-Paradigm Meta-Theorems Governing All 200 Clustering Algorithms",
        "instructions": (
            "Write a deeply rigorous 12-14 page IEEE theoretical chapter establishing the fundamental mathematical limits that constrain ALL 200 algorithms:\n"
            "1. *Kleinberg's Impossibility Theorem* (@kleinberg2002impossibility): State Scale-Invariance, Richness, and Consistency axioms formally and give a COMPLETE, self-contained proof that no clustering function satisfies all three simultaneously, followed by a table classifying which axiom each of our 8 algorithm families sacrifices!\n"
            "2. *Ackerman–Ben-David Axiomatic Relaxation & Carlsson–Mémoli Functoriality* (@ackerman2008measures, @carlsson2010characterization): Prove how relaxing clustering functions to clustering quality measures or dendrogram functors restores consistency and uniquely characterizes Single Linkage among excisive hierarchical functors.\n"
            "3. *High-Dimensional Distance Concentration (Curse of Dimensionality)* (@beyer1999nearest): Give a complete proof via Chebyshev's inequality that $(D_\"max\" - D_\"min\") / D_\"min\" ->_p 0$ as $d -> infinity$ under i.i.d. coordinates, and prove why subspace (Cat. VI), cosine/vMF (Algs. 12, 89), and SNN (Algs. 36, 52) survive while Euclidean $k$-means and fixed-$epsilon$ DBSCAN collapse.\n"
            "4. *Stability Theory & Model Selection Pitfalls* (@ben-david2006sober): Prove Ben-David et al.'s theorem that empirical clustering stability converges to symmetry of optimal minimizers rather than correctness of $k$.\n"
            "5. *Comparative Finite-Sample Breakdown Point Hierarchy* [Original Theorem & Proof]: Formally prove and rank the exact Donoho-Huber breakdown points across Centroid ($1/n$), Trimmed ($alpha$), Medoid/K-Medians ($approx 1/(2k)$), Student-t/PCM (bounded influence), and Density/Noise-isolating methods (up to $1 - \"minPts\"/n$).\n"
            "6. *Internal & External Validation Geometry* (@rousseeuw1987silhouettes, @davies1979cluster, @calinski1974dendrite, @dunn1974well, @hubert1985comparing): Prove why Silhouette/CHI/DBI are mathematically biased toward spherical convex clusters and penalize true non-convex manifold recoveries."
        ),
    },
    {
        "cid": 12, "roman": "XII", "min_p": 12, "max_p": 22,
        "title": "Complete 200-Algorithm Master Comparative Reference Matrices",
        "scope": "Exhaustive 10-Dimension Theoretical Reference Compendium for Algorithms 1-200",
        "instructions": (
            "You MUST generate `chapters/cmp_12.typ` containing the COMPLETE, exhaustive, 200-row Master Comparative Reference Tables covering EVERY SINGLE ALGORITHM from ID 1 to ID 200 (zero gaps, organized into 8 clean IEEE booktabs tables—one per category—plus a cross-category executive summary table and complexity-class taxonomy!).\n"
            "You may write a tiny helper script or write directly, ensuring EVERY SINGLE ONE of the 200 algorithms in `C:\\dev\\clusterOverview\\monograph\\data\\catalog_200.json` appears with accurate, algorithm-specific entries across columns: "
            "`ID`, `Algorithm`, `Time Complexity`, `Geometry / Bias`, `Noise / Outlier`, `Auto-k`, `High-d`, `Primary Domain`.\n"
            "Use `#table` with `text(size: 6.8pt)` and `breakable: true` (inside single-column `set page(columns: 1)` or clean full-width spanning tables, wait: NEVER call `set page` mid-chapter unless you restore `columns: 2` cleanly; even better: use `#colbreak()` / compact 6-column full-column tables OR `place(top, scope: \"parent\")`? Wait: `place` cannot break across pages! In Typst, to have multi-page full-width wide tables cleanly, simply wrap the table section in `#columns(1)[ ... ]`? NO, inside a 2-column page, `set page(columns: 1)` starts a new full-width page and `set page(columns: 2)` at the end of `cmp_12.typ` restores two-column IEEE mode seamlessly! TEST THIS with `verify_chapter.py` to be 100% sure!).\n"
            "Also include a deep 2-page mathematical commentary analyzing complexity-class stratification ($O(n)$, $O(n log n)$, $O(n k d T)$, $O(n^2)$, $O(n^3)$), deterministic vs. stochastic trade-offs, and memory footprints across all 200 methods."
        ),
    },
    {
        "cid": 13, "roman": "XIII", "min_p": 12, "max_p": 20,
        "title": "Empirical Multi-Geometry Benchmark Comparison of All 200 Algorithms",
        "scope": "Real Experimental Evaluation Across 4 Benchmark Datasets via the Native Rust ClusterOverview Engine",
        "instructions": (
            "Read `C:\\dev\\clusterOverview\\monograph\\data\\benchmark_summary.json` (and optionally the `bench_*.csv` files in `monograph/data/`) and author a deeply analytical, 12-16 page IEEE empirical evaluation chapter grounded 100% in our REAL measured results across all 4 × 200 = 800 runs!\n"
            "MUST include:\n"
            "1. Experimental protocol description (4 datasets: Iris convex/overlapping $n=120, d=4$; Customers Mixed numeric+categorical $n=150, d=6$; Multi-Shape non-convex moons+blob+ellipse+outliers $n=188, d=2$; High-D Subspace+8 noise dims+outliers $n=180, d=12$; standardized Z-score + deterministic seed 42 in native Rust `ClusterOverview`).\n"
            "2. Embed all 3 REAL vector SVG figures generated from our benchmark runs:\n"
            "   - `#wide_figure(caption: [Mean Silhouette score by category across four benchmark geometries (N = 200 algorithms).])[#image(\"../figures/fig_cat_silhouette.svg\", width: 96%)]`\n"
            "   - `#wide_figure(caption: [Quality-speed Pareto landscape across all 200 clustering algorithms (mean Silhouette vs. log10 runtime).])[#image(\"../figures/fig_ari_vs_runtime.svg\", width: 96%)]`\n"
            "   - `#wide_figure(caption: [Category Mean and Max Adjusted Rand Index (ARI) vs. ground truth across diverse geometries.])[#image(\"../figures/fig_geometry_heatmap.svg\", width: 96%)]`\n"
            "3. Complete Category-Level Empirical Performance Table (exact real Mean/Max Silhouette, Mean/Max ARI, Mean/Median ms for all 8 categories across all 4 datasets).\n"
            "4. Top-25 Overall & Per-Geometry Leaderboard Tables with exact real numbers.\n"
            "5. Complete 200-Algorithm Empirical Benchmark Appendix Table (all 200 IDs with real Iris Sil/ARI, Customers Sil, Multi-Shape Sil/ARI, High-D Sil/ARI, Mean ms — cleanly formatted using `set page(columns: 1)` and restored with `set page(columns: 2)` at the end!).\n"
            "6. Rigorous theoretical-empirical reconciliation: formally connect the observed performance shifts (e.g., why Silhouette and ARI diverge on non-convex Multi-Shape, why distance concentration hits Euclidean methods on High-D, why runtime spans 4 orders of magnitude from $<0.1$ ms to $>20$ ms) directly back to our theorems in Parts I-VIII and Comparison Chapter XI!"
        ),
    },
    {
        "cid": 14, "roman": "XIV", "min_p": 10, "max_p": 16,
        "title": "Practitioner's Engineering Applications Handbook and Algorithm Selection Architecture",
        "scope": "Domain-to-Algorithm Mathematical Selection Framework Across 10 Scientific & Industrial Verticals",
        "instructions": (
            "Write a comprehensive 12-14 page IEEE applications handbook synthesizing the engineering and scientific deployment of all 200 clustering algorithms across 10 major verticals:\n"
            "1. *Genomics, Transcriptomics & Systems Biology* (scRNA-seq, gene expression biclustering, phylogenetics: Leiden, Louvain, UPGMA/NJ, Cheng-Church, FABIA, Consensus Monti, HDBSCAN)\n"
            "2. *Computer Vision, Medical Imaging & 3D LiDAR Robotics* (SLIC, Watershed, Mean Shift++, Ncut, DBSCAN-3D, sFCM, SwAV, SCAN, Contrastive Clustering)\n"
            "3. *Natural Language Processing, Retrieval & Large-Scale Embeddings* (Spherical K-Means, vMF, LDA, HDP, CTM, DMM, BTM, K-Means++ IVF-PQ, Chinese Whispers)\n"
            "4. *Cybersecurity, Fraud & Anomaly Detection* (OPTICS-OF, DBSCAN, HDBSCAN, DenStream, CluStream, D-Stream, Trimmed K-Means, Student-t TMM)\n"
            "5. *Quantitative Finance, Portfolio Risk & Customer Analytics* (K-Prototypes, PAM/FasterPAM, MCLUST, Ward, W-K-Means, k-Shape, Soft-DTW)\n"
            "6. *High-Velocity IoT, Network Telemetry & Edge Streaming* (BIRCH, BICO, StreamKM++, ClusTree, CluStream, DenStream, D-Stream, Mini-Batch K-Means)\n"
            "7. *Geospatial Analytics, Urban Mobility & Spatiotemporal Tracking* (ST-DBSCAN, TRACLUS, CLARANS, STING, WaveCluster, OPTICS, VDBSCAN)\n"
            "8. *Social, Web, Communication & Biological Network Science* (Leiden, Louvain, Infomap, MCL, SBM/DC-SBM, BigCLAM, OSLOM, CPM, METIS)\n"
            "9. *High-Dimensional Industrial Sensors & Subspace Discovery* (SSC, PROCLUS, ORCLUS, SUBCLU, CLIQUE, COSA, P3C, STATPC)\n"
            "10. *Multi-Modal Fusion & Unsupervised Consensus Pipelines* (Multi-View K-Means, Co-regularized Spectral, CSPA, HGPA, MCLA, VaDE, IDEC).\n"
            "Include formal *Algorithm Selection Decision Theorems / Rules* (mapping observable data invariants: $n$, $d$, metric/non-metric, sparsity, outlier fraction, cluster geometry, streaming budget, interpretability -> provably appropriate candidate families), "
            "a Master Diagnostic Decision Matrix Table, and a Production Engineering Checklist (preprocessing, scaling, validation pitfalls, and multi-algorithm ensemble verification using ClusterOverview)."
        ),
    },
]


def build_cmp_prompt(spec: dict) -> str:
    cid = spec["cid"]
    target = f"C:\\dev\\clusterOverview\\monograph\\chapters\\cmp_{cid:02d}.typ"
    return f"""You are a senior applied mathematician and IEEE Transactions lead author contributing Comparison Chapter {spec['roman']} to Part IX of the monograph "Mathematical Foundations, Formal Proofs, and Applications of 200 Clustering Algorithms" (Created by Gemini, AI made by Google).

TARGET FILE: `{target}`
PAGE TARGET: **{spec['min_p']} to {spec['max_p']} pages** of rich, authoritative, zero-filler IEEE two-column mathematical prose, theorems, step-by-step proofs, and tables!

HEADER FORMAT (MUST start file with):
```typst
#import "../template/ieee.typ": *

#cmp_entry(
  roman: "{spec['roman']}",
  title: "{spec['title']}",
  scope_label: "{spec['scope']}",
  provenance: "Original Unifying Proofs + Cited Classical Theorems",
)
```

DETAILED CHAPTER REQUIREMENTS:
{spec['instructions']}

STRICT TYPST 0.15 SYNTAX & ANTI-OVERFLOW RULES:
- Bold is SINGLE `*bold*` (NEVER `**bold**`!). Italic is `_italic_`.
- Inner products: `chevron.l u, v chevron.r` (NEVER `angle.l`/`angle.r`!).
- Subscripts before parens MUST have a space: `D_phi (x || y)`, `hat(f)_h (x)`, `lambda_k (L)`.
- Literal hash sign: write `Alg. 12` or `\\#12` (NEVER bare `#12`!).
- Multi-letter math identifiers/operators MUST be quoted: `"Tr"`, `"Var"`, `"Cov"`, `"diag"`, `"KL"`, `"cut"`, `"vol"`, `"deg"`, `"rank"`, `"span"`, `"dist"`, `"med"`, `"sign"`, `"poly"`, `"supp"`, `"ELBO"`, `"BIC"`, `"MDL"`, `"ARI"`, `"SSE"`, `"SBD"`, `"DTW"`, `"Ncut"`, `"Rcut"`, `"kNN"`, `"SNN"`, `"Gini"`, `"LOF"`. Use `arg min` and `arg max`.
- Set, order & calculus symbols in Typst 0.15: intersection is `inter` (NEVER `sect`!), big product is `product` (NEVER `prod`!), PSD order is `succ.eq` / `prec.eq` (NEVER `succeq`/`preceq`!), equivalence is `tilde.op` (NEVER `sim`!), direct/Minkowski sum is `plus.o` (NEVER `plus.circle`!), membership negation is `in.not` (NEVER `notin`!), differential is `dif` / `partial` (NEVER `diff`!). NEVER use LaTeX backslash macros like `\\text{...}`, `\\epsilon`, `\\lambda`, `\\ell`, `\\omega`. Escape literal asterisks (`DBSCAN\\*`, `A\\*`) and put spaces around `/` next to bold/italic (`*a* / *b*`, NEVER `*a*/*b*`).
- Split long display equations across 2-3 lines (`&` and `\\`) so NOTHING overflows a 3.5-inch IEEE column!
- Keep `origin:` badge text short (`"Original Proof"`, `"Original Unified Proof"`, `"Original Counterexample"`, `"Classical — Cited"`).
- If you temporarily switch to `#set page(columns: 1)` for a multi-page wide master table (Chapters XII or XIII), you MUST restore `#set page(columns: 2)` at the very end of the file!
- ONLY cite verified BibTeX keys from this exact list (NEVER invent undefined `@keys`!):
[{BIB_KEYS_STR}]

VERIFICATION:
After writing `{target}`, run:
`python C:\\dev\\clusterOverview\\monograph\\scripts\\verify_chapter.py chapters\\cmp_{cid:02d}.typ {spec['min_p']} {spec['max_p']}`
Fix ANY compile errors, warnings, or page-count shortfalls (`content_pages < {spec['min_p']}`) until it returns `"ok": true`, `"length_ok": true`, `"warnings": []`!
Then call `finish` with the verified page count.
"""


CMP_SCHEMA = {
    "type": "object",
    "properties": {
        "cid": {"type": "integer"},
        "pages": {"type": "integer"},
        "summary": {"type": "string"},
    },
    "required": ["cid", "pages"],
}


quota_hit = False


async def run_cmp_safe(spec: dict, attempt: int = 1):
    global quota_hit
    cid, lo, hi = spec["cid"], spec["min_p"], spec["max_p"]
    if quota_hit:
        return {"cid": cid, "pages": 0, "ok": False}
    ok, pgs = is_cmp_verified(cid, lo, hi)
    if ok:
        log(f"SKIP cmp_{cid:02d} — already verified ({pgs} pages)")
        return {"cid": cid, "pages": pgs, "ok": True}

    role = f"Comparison Ch. {spec['roman']}" + (f" (retry {attempt})" if attempt > 1 else "")
    try:
        await agent(
            build_cmp_prompt(spec),
            schema=CMP_SCHEMA,
            role=role,
            type_name="self",
            workspace="inherit",
        )
        ok2, pgs2 = is_cmp_verified(cid, lo, hi)
        if not ok2 and attempt < 2 and not quota_hit:
            log(f"RETRY cmp_{cid:02d} (got ok={ok2}, pages={pgs2})")
            return await run_cmp_safe(spec, attempt=attempt + 1)
        log(f"CMP_{cid:02d} -> ok={ok2}, pages={pgs2}")
        return {"cid": cid, "pages": pgs2, "ok": ok2}
    except Exception as exc:
        msg = str(exc)
        if "Individual quota reached" in msg:
            quota_hit = True
            log(f"INDIVIDUAL QUOTA HIT cmp_{cid:02d} — halting cleanly!")
        elif ("429" in msg or "ResourceExhausted" in msg) and attempt < 3:
            wait_s = 20 * attempt
            log(f"TRANSIENT 429 cmp_{cid:02d} (attempt {attempt}) — sleeping {wait_s}s and retrying...")
            await asyncio.sleep(wait_s)
            return await run_cmp_safe(spec, attempt=attempt + 1)
        elif "429" in msg or "ResourceExhausted" in msg:
            quota_hit = True
            log(f"PERSISTENT 429 cmp_{cid:02d} — halting cleanly!")
        else:
            log(f"ERROR cmp_{cid:02d} (attempt {attempt}): {msg[:300]}")
            if attempt < 2 and not quota_hit:
                return await run_cmp_safe(spec, attempt=attempt + 1)
        ok_post, pgs_post = is_cmp_verified(cid, lo, hi)
        return {"cid": cid, "pages": pgs_post, "ok": ok_post}


phase("Part IX: 14 Cross-Algorithm Comparison & Unification Chapters")
results = []
for w in range(0, len(CMP_SPECS), wave_size):
    if quota_hit:
        for rem in CMP_SPECS[w:]:
            results.append({"cid": rem["cid"], "pages": 0, "ok": False})
        break
    wave = CMP_SPECS[w : w + wave_size]
    wres = await parallel(wave, lambda s: run_cmp_safe(s))
    results.extend(wres)

ok_cnt = sum(1 for r in results if r["ok"])
tot_pages = sum(r["pages"] for r in results)
failed = [r["cid"] for r in results if not r["ok"]]
log(f"COMPARISONS DONE: {ok_cnt}/{len(CMP_SPECS)} verified, failed={failed}, total_pages={tot_pages}")
print(json.dumps({"verified": ok_cnt, "failed": failed, "total_pages": tot_pages}, indent=2))
