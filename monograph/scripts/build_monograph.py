import json
import subprocess
import sys
from pathlib import Path

MONO = Path(r"C:\dev\clusterOverview\monograph")
CH_DIR = MONO / "chapters"
catalog = json.loads((MONO / "data" / "catalog_200.json").read_text(encoding="utf-8"))

PARTS = [
    (
        1, "I", "Centroid, Medoid, and Partitioning Algorithms", "Algorithms 1–30",
        "Part I examines thirty prototype-based partitioning algorithms spanning Euclidean, Manhattan, spherical, Bregman, kernelized, trimmed, weighted, constrained, deterministic, and parallelizable vector-quantization regimes.",
        "Every method in Part I minimizes an empirical distortion functional over a discrete Voronoi or exemplar partition, unified by alternating projection onto prototype manifolds and convex or Bregman cell decompositions.",
    ),
    (
        2, "II", "Density-Based and Mode-Seeking Algorithms", "Algorithms 31–55",
        "Part II investigates twenty-five non-parametric density-connectivity, level-set, kernel gradient-ascent, and density-peak algorithms capable of recovering arbitrary non-convex manifolds while isolating background noise.",
        "Clusters are formulated either as connected components of upper density level sets {x : f(x) >= lambda} or as basins of attraction of local maxima under continuous or discrete kernel density gradient flows.",
    ),
    (
        3, "III", "Hierarchical and Tree-Based Algorithms", "Algorithms 56–80",
        "Part III analyzes twenty-five agglomerative, divisive, metric-tree, phylogenetic, and extreme-scale hierarchical clustering methods constructing multi-resolution nested taxonomies.",
        "Hierarchical clustering maps finite dissimilarity spaces into ultrametric trees (dendrograms) governed by the Lance–Williams recursion, minimum spanning tree equivalence, or hierarchical coreset summarization.",
    ),
    (
        4, "IV", "Distribution, Statistical Mixture, and Topic Models", "Algorithms 81–102",
        "Part IV develops twenty-two generative probabilistic mixture models, Bayesian nonparametric processes, conceptual utility trees, and latent topic/factorization models.",
        "Observations arise from latent-variable mixture densities maximized via Expectation–Maximization, variational lower bounds (ELBO), collapsed Gibbs sampling, Bayesian Information Criterion, or Minimum Message Length.",
    ),
    (
        5, "V", "Graph, Spectral, and Network Community Detection", "Algorithms 103–132",
        "Part V covers thirty spectral graph cut, modularity optimization, random-walk flow, message-passing, overlapping community, multilevel partitioning, and stochastic blockmodel algorithms.",
        "Partitions optimize discrete graph cut ratios, Newman–Girvan modularity, Markov random-walk description length, or generative block-edge likelihoods via Laplacian eigen-relaxations and local greedy/message updates.",
    ),
    (
        6, "VI", "Grid, Subspace, and Projected Clustering Algorithms", "Algorithms 133–150",
        "Part VI treats eighteen spatial grid quantization, wavelet multiresolution, axis-parallel subspace, arbitrarily oriented projected, and sparse self-expressive subspace clustering algorithms designed to defeat the curse of dimensionality.",
        "High-dimensional structure is recovered by restricting distortion or density evaluation to cluster-specific low-dimensional subspaces via downward-closed pruning, entropy criteria, local PCA, or L1 sparse representation.",
    ),
    (
        7, "VII", "Fuzzy, Biclustering, Time-Series, and Domain-Specific Methods", "Algorithms 151–175",
        "Part VII addresses twenty-five soft-membership fuzzy/possibilistic partitioners, two-mode biclustering models, elastic/shape-invariant time-series clusterers, superpixel/3D/trajectory algorithms, consensus ensembles, and multi-view methods.",
        "Domain invariances—gradual membership, simultaneous row-column coherence, phase/warp invariance, spatial continuity, and cross-view consensus—are encoded directly into generalized dissimilarity functionals and hypergraph objectives.",
    ),
    (
        8, "VIII", "Neural, Deep, Streaming, and Evolutionary Algorithms", "Algorithms 176–200",
        "Part VIII studies twenty-five competitive neural maps, deep representation-clustering architectures, single-pass evolving data-stream engines, and population-based stochastic global metaheuristics.",
        "Clustering is coupled either with nonlinear representation learning and self-supervised mutual-information/optimal-transport losses, bounded-memory fading micro-cluster coresets, or global stochastic search dynamics.",
    ),
]


def esc(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("#", "\\#")


def build_master_typ():
    missing_algs = [i for i in range(1, 201) if not (CH_DIR / f"alg_{i:03d}.typ").exists()]
    missing_cmps = [i for i in range(1, 15) if not (CH_DIR / f"cmp_{i:02d}.typ").exists()]
    if missing_algs or missing_cmps:
        print(f"WARNING: missing algs={missing_algs}, missing cmps={missing_cmps}")

    lines = [
        '#import "template/ieee.typ": *',
        '#show: ieee_doc.with()',
        '',
        '// Frontmatter: IEEE Title, Abstract, Notation, Provenance & TOC',
        '#include "frontmatter.typ"',
        '',
    ]

    for cnum, roman, title, arange, summary, principle in PARTS:
        cat_items = [it for it in catalog if it["cat_num"] == cnum]
        # Build crisp 2-column grid listing all algorithms in this Part
        half = (len(cat_items) + 1) // 2
        left_col = " \\ ".join(f"*\\#{it['id']}.* {esc(it['name'])}" for it in cat_items[:half])
        right_col = " \\ ".join(f"*\\#{it['id']}.* {esc(it['name'])}" for it in cat_items[half:])
        lines.append(f'#part_opener(')
        lines.append(f'  roman: "{roman}",')
        lines.append(f'  title: "{esc(title)}",')
        lines.append(f'  alg_range: "{arange}",')
        lines.append(f'  summary: [{esc(summary)}],')
        lines.append(f'  core_principle: [{esc(principle)}],')
        lines.append(f'  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [{left_col}], [{right_col}]),')
        lines.append(f')')
        for it in cat_items:
            p = CH_DIR / f"alg_{it['id']:03d}.typ"
            if p.exists():
                lines.append(f'#include "chapters/alg_{it["id"]:03d}.typ"')
        lines.append('')

    # Part IX: Comparisons
    cmp_titles = [
        "I. Grand Unified Variational, Geometric & Information-Theoretic Framework",
        "II. Formal Mathematical Equivalences, Asymptotic Limits & Reductions",
        "III. Intra-Category Comparison: Centroid, Medoid & Partitioning (1–30)",
        "IV. Intra-Category Comparison: Density & Mode-Seeking (31–55)",
        "V. Intra-Category Comparison: Hierarchical & Tree-Based (56–80)",
        "VI. Intra-Category Comparison: Distribution, Statistical & Topic Models (81–102)",
        "VII. Intra-Category Comparison: Graph, Spectral & Network (103–132)",
        "VIII. Intra-Category Comparison: Grid, Subspace & Projection (133–150)",
        "IX. Intra-Category Comparison: Fuzzy, Biclustering, Time-Series & Domain (151–175)",
        "X. Intra-Category Comparison: Neural, Deep, Streaming & Evolutionary (176–200)",
        "XI. Axiomatic Foundations, Impossibility Theorems & Breakdown Robustness",
        "XII. Complete 200-Algorithm Master Comparative Reference Matrices",
        "XIII. Empirical Multi-Geometry Benchmark Comparison (ClusterOverview Engine)",
        "XIV. Practitioner's Engineering Applications Handbook & Selection Architecture",
    ]
    left_c = " \\ ".join(esc(t) for t in cmp_titles[:7])
    right_c = " \\ ".join(esc(t) for t in cmp_titles[7:])
    lines.extend([
        '#part_opener(',
        '  roman: "IX",',
        '  title: "Cross-Algorithm Mathematical Unification, Equivalences, and Empirical Comparisons",',
        '  alg_range: "Comparative Chapters I–XIV",',
        '  summary: [Part IX synthesizes the entire 200-algorithm compendium through fourteen dedicated comparative chapters proving cross-family mathematical equivalences, axiomatic impossibility boundaries, intra-category head-to-head analyses, master theoretical reference matrices, real multi-geometry benchmarks, and domain engineering decision architectures.],',
        '  core_principle: [Apparent algorithmic diversity collapses into a tightly connected web of variational limits, spectral-kernel dualities, ultrametric level-set functors, and information-theoretic trade-offs governed by Kleinberg impossibility and high-dimensional measure concentration.],',
        f'  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [{left_c}], [{right_c}]),',
        ')',
    ])
    for cid in range(1, 15):
        p = CH_DIR / f"cmp_{cid:02d}.typ"
        if p.exists():
            lines.append(f'#include "chapters/cmp_{cid:02d}.typ"')

    lines.extend([
        '',
        '#pagebreak(weak: true)',
        '#place(top + center, scope: "parent", float: true, block(width: 100%, fill: ieee_navy, inset: (x: 16pt, y: 12pt), radius: 3pt, below: 10pt, {',
        '  set align(center)',
        '  heading(level: 1, outlined: true, numbering: none)[Complete IEEE Bibliography]',
        '}))',
        '#bibliography("references.bib", title: none, style: "ieee")',
    ])

    main_typ = MONO / "main.typ"
    main_typ.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {main_typ} (algs present: {200-len(missing_algs)}/200, cmps present: {14-len(missing_cmps)}/14)")
    return main_typ


def compile_pdf(main_typ: Path):
    out_pdf = MONO / "ClusterOverview_200_Clustering_Algorithms_IEEE_Monograph.pdf"
    res = subprocess.run(
        ["typst", "compile", str(main_typ), str(out_pdf)],
        cwd=str(MONO),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if res.returncode != 0:
        print("COMPILE FAILED:\n", res.stderr[:6000])
        sys.exit(1)
    if res.stderr.strip():
        print("WARNINGS:\n", res.stderr[:3000])
    # Count exact pages from PDF
    import re
    pdf_bytes = out_pdf.read_bytes()
    matches = re.findall(rb"/Type\s*/Pages\b[^>]*?/Count\s+(\d+)", pdf_bytes, flags=re.DOTALL)
    pages = max(int(m) for m in matches) if matches else 0
    size_mb = round(len(pdf_bytes) / (1024 * 1024), 2)
    print(f"SUCCESS! Compiled {out_pdf} -> {pages} pages ({size_mb} MB)")


if __name__ == "__main__":
    mt = build_master_typ()
    if "--compile" in sys.argv:
        compile_pdf(mt)
