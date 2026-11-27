import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from pypdf import PdfReader, PdfWriter

MONO = Path(__file__).resolve().parent.parent
CH_DIR = MONO / "chapters"
TMP_DIR = MONO / "_vol_build"
TMP_DIR.mkdir(parents=True, exist_ok=True)
QA_DIR = MONO / "_qa_preview"
QA_DIR.mkdir(parents=True, exist_ok=True)

catalog = json.loads((MONO / "data" / "catalog_200.json").read_text(encoding="utf-8"))
bib_keys = json.loads((MONO / "data" / "bib_keys.json").read_text(encoding="utf-8"))

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

CMP_TITLES = [
    (1, "I", "Grand Unified Variational, Geometric & Information-Theoretic Framework"),
    (2, "II", "Formal Mathematical Equivalences, Asymptotic Limits & Reductions"),
    (3, "III", "Intra-Category Comparison: Centroid, Medoid & Partitioning (1–30)"),
    (4, "IV", "Intra-Category Comparison: Density & Mode-Seeking (31–55)"),
    (5, "V", "Intra-Category Comparison: Hierarchical & Tree-Based (56–80)"),
    (6, "VI", "Intra-Category Comparison: Distribution, Statistical & Topic Models (81–102)"),
    (7, "VII", "Intra-Category Comparison: Graph, Spectral & Network (103–132)"),
    (8, "VIII", "Intra-Category Comparison: Grid, Subspace & Projection (133–150)"),
    (9, "IX", "Intra-Category Comparison: Fuzzy, Biclustering, Time-Series & Domain (151–175)"),
    (10, "X", "Intra-Category Comparison: Neural, Deep, Streaming & Evolutionary (176–200)"),
    (11, "XI", "Axiomatic Foundations, Impossibility Theorems & Breakdown Robustness"),
    (12, "XII", "Complete 200-Algorithm Master Comparative Reference Matrices"),
    (13, "XIII", "Empirical Multi-Geometry Benchmark Comparison (ClusterOverview Engine)"),
    (14, "XIV", "Practitioner's Engineering Applications Handbook & Selection Architecture"),
]

# Deterministic canonical order for all 244 verified BibTeX keys so [1]..[244] match identically in ALL volumes!
HIDDEN_CITE_BLOCK = (
    '#place(top + left, float: false, hide(block(width: 0pt, height: 0pt, clip: true, { set text(size: 0.1pt); ['
    + " ".join(f"@{k}" for k in bib_keys)
    + '] })))'
)


def esc(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("#", "\\#")


def make_part_opener(cnum: int) -> str:
    _, roman, title, arange, summary, principle = PARTS[cnum - 1]
    cat_items = [it for it in catalog if it["cat_num"] == cnum]
    half = (len(cat_items) + 1) // 2
    left_col = " \\ ".join(f"*\\#{it['id']}.* {esc(it['name'])}" for it in cat_items[:half])
    right_col = " \\ ".join(f"*\\#{it['id']}.* {esc(it['name'])}" for it in cat_items[half:])
    return "\n".join([
        "#part_opener(",
        f'  roman: "{roman}",',
        f'  title: "{esc(title)}",',
        f'  alg_range: "{arange}",',
        f'  summary: [{esc(summary)}],',
        f'  core_principle: [{esc(principle)}],',
        f'  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [{left_col}], [{right_col}]),',
        ")",
    ])


def make_part9_opener() -> str:
    left_c = " \\ ".join(esc(f"{r}. {t}") for _, r, t in CMP_TITLES[:7])
    right_c = " \\ ".join(esc(f"{r}. {t}") for _, r, t in CMP_TITLES[7:])
    return "\n".join([
        "#part_opener(",
        '  roman: "IX",',
        '  title: "Cross-Algorithm Mathematical Unification, Equivalences, and Empirical Comparisons",',
        '  alg_range: "Comparative Chapters I–XIV",',
        '  summary: [Part IX synthesizes the entire 200-algorithm compendium through fourteen dedicated comparative chapters proving cross-family mathematical equivalences, axiomatic impossibility boundaries, intra-category head-to-head analyses, master theoretical reference matrices, real multi-geometry benchmarks, and domain engineering decision architectures.],',
        '  core_principle: [Apparent algorithmic diversity collapses into a tightly connected web of variational limits, spectral-kernel dualities, ultrametric level-set functors, and information-theoretic trade-offs governed by Kleinberg impossibility and high-dimensional measure concentration.],',
        f'  alg_list: grid(columns: (1fr, 1fr), column-gutter: 12pt, [{left_c}], [{right_c}]),',
        '  list_heading: "COMPARATIVE CHAPTERS IN PART IX",',
        ")",
    ])


def build_volume_typ(vol_name: str, body_lines: list, start_page: int, show_bib: bool) -> Path:
    # Put HIDDEN_CITE_BLOCK inside page background/after first opener? Even better:
    # placing it inside the VERY FIRST page header or right AFTER body_lines keeps citation order?
    # Wait: IEEE style numbers citations in order of FIRST appearance! So to keep [1]..[244] identical
    # WITHOUT breaking #part_opener's leading #pagebreak(weak: true), we put HIDDEN_CITE_BLOCK
    # inside a tiny show-rule or simply omit #pagebreak(weak: true) on the VERY FIRST opener!
    # Even simpler: emit #counter(page).update(start_page) + body_lines (first line is #part_opener),
    # and place HIDDEN_CITE_BLOCK right AFTER the first #part_opener (before any chapter cites anything!).
    lines = [
        '#import "../template/ieee.typ": *',
        '#show: ieee_doc.with()',
        f'#counter(page).update({start_page})',
        "",
    ]
    # Insert HIDDEN_CITE_BLOCK immediately after the first #part_opener(...) block (before any @citation)
    if body_lines:
        lines.append(body_lines[0])
        lines.append(HIDDEN_CITE_BLOCK)
        lines.extend(body_lines[1:])
    lines.append("")
    if show_bib:
        lines.extend([
            '#pagebreak(weak: true)',
            '#place(top + center, scope: "parent", float: true, block(width: 100%, fill: ieee_navy, inset: (x: 16pt, y: 12pt), radius: 3pt, below: 10pt, {',
            '  set align(center)',
            '  heading(level: 1, outlined: true, numbering: none)[Complete IEEE Bibliography]',
            '}))',
            '#context [#metadata((kind: "part", label: "Complete IEEE Bibliography", page: counter(page).get().first())) <toc_mark>]',
            '#bibliography("../references.bib", title: none, full: true, style: "ieee")',
        ])
    else:
        lines.extend([
            '#show bibliography: none',
            '#bibliography("../references.bib", full: true, style: "ieee")',
        ])
    p = TMP_DIR / f"{vol_name}.typ"
    p.write_text("\n".join(lines), encoding="utf-8")
    return p


def compile_and_query_vol(typ_path: Path, pdf_path: Path):
    cmd = ["typst", "compile", "--root", str(MONO), "--jobs", "2", str(typ_path), str(pdf_path)]
    res = subprocess.run(cmd, cwd=str(MONO), capture_output=True, text=True, encoding="utf-8", errors="replace")
    if res.returncode != 0:
        print(f"FAILED compiling {typ_path.name}:\n{res.stderr[:4000]}")
        sys.exit(1)
    warns = [ln for ln in res.stderr.splitlines() if ln.strip()]
    if warns:
        print(f"  [warnings in {typ_path.name}]: {warns[:5]}")
    reader = PdfReader(str(pdf_path))
    n_pages = len(reader.pages)
    q = subprocess.run(
        ["typst", "query", "--root", str(MONO), str(typ_path), "<toc_mark>", "--field", "value"],
        cwd=str(MONO),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    marks = json.loads(q.stdout) if (q.returncode == 0 and q.stdout.strip()) else []
    return n_pages, marks


def make_vol_bodies():
    # Split into 5 balanced volumes (~195-225 pages each, <0.9 GB RAM each!)
    vol_defs = [
        ("vol1_parts1_2", [1, 2], False),
        ("vol2_parts3_4", [3, 4], False),
        ("vol3_parts5_6", [5, 6], False),
        ("vol4_parts7_8", [7, 8], False),
        ("vol5_part9_bib", [9], True),
    ]
    bodies = []
    for vname, parts_list, show_bib in vol_defs:
        lines = []
        for pnum in parts_list:
            if pnum <= 8:
                _, roman, title, _, _, _ = PARTS[pnum - 1]
                lines.append(make_part_opener(pnum))
                lines.append(f'#context [#metadata((kind: "part", roman: "{roman}", label: "Part {roman}: {esc(title)}", page: counter(page).get().first())) <toc_mark>]')
                for it in [x for x in catalog if x["cat_num"] == pnum]:
                    lines.append(f'#include "../chapters/alg_{it["id"]:03d}.typ"')
                    lines.append(f'#context [#metadata((kind: "alg", part: "{roman}", id: {it["id"]}, label: "Algorithm {it["id"]}: {esc(it["name"])}", page: counter(page).get().first())) <toc_mark>]')
            else:
                lines.append(make_part9_opener())
                lines.append('#context [#metadata((kind: "part", roman: "IX", label: "Part IX: Cross-Algorithm Mathematical Unification & Comparisons", page: counter(page).get().first())) <toc_mark>]')
                for cid, roman, title in CMP_TITLES:
                    lines.append(f'#include "../chapters/cmp_{cid:02d}.typ"')
                    lines.append(f'#context [#metadata((kind: "cmp", part: "IX", id: {cid}, label: "Comparison Ch. {roman}: {esc(title)}", page: counter(page).get().first())) <toc_mark>]')
        bodies.append((vname, lines, show_bib))
    return bodies


def extract_entry_start_pages(typ_path: Path, start_page: int):
    """Query exact heading start pages (level 1 and level 2) via Typst metadata placed BEFORE each chapter?
    Actually, querying `heading.where(level: 1)` / `level: 2` locations gives EXACT start pages!"""
    # Query exact start page of every level-1 and level-2 heading by appending a context probe?
    # Even simpler & exact: typst query outputs location page numbers when we emit metadata with `here()`!
    pass


def make_vol_bodies_exact():
    """Place `<toc_mark>` AFTER the pagebreak/banner by querying heading locations cleanly inside Typst!"""
    vol_defs = [
        ("vol1_parts1_2", [1, 2], False),
        ("vol2_parts3_4", [3, 4], False),
        ("vol3_parts5_6", [5, 6], False),
        ("vol4_parts7_8", [7, 8], False),
        ("vol5_part9_bib", [9], True),
    ]
    bodies = []
    for vname, parts_list, show_bib in vol_defs:
        lines = []
        for pnum in parts_list:
            if pnum <= 8:
                lines.append(make_part_opener(pnum))
                for it in [x for x in catalog if x["cat_num"] == pnum]:
                    lines.append(f'#include "../chapters/alg_{it["id"]:03d}.typ"')
            else:
                lines.append(make_part9_opener())
                for cid, _, _ in CMP_TITLES:
                    lines.append(f'#include "../chapters/cmp_{cid:02d}.typ"')
        # Emit ONE single metadata array at end of volume querying EVERY level 1 & level 2 heading's exact page!
        lines.append('''#context [
  #let hs = query(selector(heading.where(level: 1)).or(heading.where(level: 2)))
  #let rows = hs.map(h => (
    level: h.level,
    title: if h.body.has("text") { h.body.text } else { repr(h.body) },
    page: counter(page).at(h.location()).first(),
  ))
  #metadata(rows) <toc_mark>
]''')
        bodies.append((vname, lines, show_bib))
    return bodies


def clean_heading_text(raw: str) -> str:
    # Typst repr(sequence(...)) fallback parser if heading body had interpolated elements
    if not raw.startswith("sequence(") and not raw.startswith("["):
        return raw.strip()
    parts = re.findall(r'\[([^\[\]]*)\]', raw)
    joined = "".join(parts).strip()
    return joined if joined else raw


def build_frontmatter_typ(toc_entries: list) -> Path:
    # Build authentic frontmatter (Title + Abstract + Provenance + Notation) + Exact Global Master TOC!
    fm_src = (MONO / "frontmatter.typ").read_text(encoding="utf-8")
    # Strip the auto #outline(...) at the bottom and replace with our exact synchronized global TOC
    cut_idx = fm_src.find("#outline(")
    if cut_idx != -1:
        fm_src = fm_src[:cut_idx]

    toc_lines = [
        '#set par(first-line-indent: 0pt, justify: false, leading: 0.48em, spacing: 0.28em)',
        '#set text(size: 7.7pt)',
    ]
    for entry in toc_entries:
        lvl = entry["level"]
        title = esc(entry["title"])
        pg = entry["page"]
        if lvl == 1:
            toc_lines.append(
                f'#v(4pt) #grid(columns: (1fr, auto), column-gutter: 6pt, text(size: 8.4pt, weight: "bold", fill: ieee_navy)[{title}], text(size: 8.4pt, weight: "bold", fill: ieee_navy)[{pg}]) #v(-2pt) #line(length: 100%, stroke: 0.35pt + ieee_border)'
            )
        else:
            toc_lines.append(
                f'#grid(columns: (1fr, auto), column-gutter: 4pt, [#h(6pt){title} #box(width: 1fr, repeat([.]))], [{pg}])'
            )

    full_fm = "\n".join([
        '#import "../template/ieee.typ": *',
        '#show: ieee_doc.with()',
        HIDDEN_CITE_BLOCK,
        fm_src.replace('#import "template/ieee.typ": *', ''),
        "\n".join(toc_lines),
        '#show bibliography: none',
        '#bibliography("../references.bib", full: true, style: "ieee")',
    ])
    p = TMP_DIR / "vol0_frontmatter.typ"
    p.write_text(full_fm, encoding="utf-8")
    return p


def main():
    bodies = make_vol_bodies_exact()

    # PASS 1: Measure exact page counts of Vols 1..5 (start at 1)
    print("=== PASS 1: Measuring exact page counts across Volumes 1-5 ===")
    vol_page_counts = []
    pass1_headings = []
    for vname, blines, show_bib in bodies:
        tpath = build_volume_typ(vname, blines, start_page=1, show_bib=show_bib)
        ppath = TMP_DIR / f"{vname}.pdf"
        n_pgs, marks = compile_and_query_vol(tpath, ppath)
        hs = marks[0] if marks else []
        print(f"  {vname}: {n_pgs} pages, {len(hs)} headings")
        vol_page_counts.append(n_pgs)
        pass1_headings.append(hs)

    # Estimate Frontmatter page count via a dry compile with provisional page offsets, then stabilize (1-2 iterations)
    print("=== PASS 2: Building exact synchronized Frontmatter & Master TOC ===")
    fm_pages = 7  # initial guess
    for _ in range(3):
        toc_entries = []
        cursor = fm_pages + 1
        for n_pgs, hs in zip(vol_page_counts, pass1_headings):
            for h in hs:
                # Note: page in Pass 1 started at 1, so global page = cursor + (h["page"] - 1)
                toc_entries.append({
                    "level": h["level"],
                    "title": clean_heading_text(h["title"]),
                    "page": cursor + (int(h["page"]) - 1),
                })
            cursor += n_pgs
        fm_typ = build_frontmatter_typ(toc_entries)
        fm_pdf = TMP_DIR / "vol0_frontmatter.pdf"
        actual_fm_pages, _ = compile_and_query_vol(fm_typ, fm_pdf)
        print(f"  Frontmatter dry-run (assumed={fm_pages}) -> actual={actual_fm_pages} pages")
        if actual_fm_pages == fm_pages:
            break
        fm_pages = actual_fm_pages

    # PASS 3: Recompile Volumes 1..5 with exact global starting page numbers!
    print("=== PASS 3: Compiling Volumes 1-5 with exact global continuous page numbers ===")
    cursor = fm_pages + 1
    final_pdfs = [TMP_DIR / "vol0_frontmatter.pdf"]
    final_toc = []
    vol_ranges = []
    for idx, (vname, blines, show_bib) in enumerate(bodies):
        tpath = build_volume_typ(vname, blines, start_page=cursor, show_bib=show_bib)
        ppath = TMP_DIR / f"{vname}.pdf"
        n_pgs, marks = compile_and_query_vol(tpath, ppath)
        hs = marks[0] if marks else []
        print(f"  {vname}: global pages {cursor}..{cursor + n_pgs - 1} ({n_pgs} pages)")
        vol_ranges.append((vname, tpath, cursor, n_pgs))
        for h in hs:
            final_toc.append({
                "level": h["level"],
                "title": clean_heading_text(h["title"]),
                "page": int(h["page"]),
            })
        final_pdfs.append(ppath)
        cursor += n_pgs

    # Rebuild Frontmatter once more if any volume shifted by +/- 1 page due to header width
    if len(final_toc) == len(toc_entries) and any(a["page"] != b["page"] for a, b in zip(final_toc, toc_entries)):
        print("  Refreshing Frontmatter TOC with exact Pass-3 page numbers...")
        fm_typ = build_frontmatter_typ(final_toc)
        fm_pdf = TMP_DIR / "vol0_frontmatter.pdf"
        actual_fm2, _ = compile_and_query_vol(fm_typ, fm_pdf)
        assert actual_fm2 == fm_pages, f"Frontmatter page shift: {fm_pages} -> {actual_fm2}"

    # PASS 4: Losslessly merge all 6 PDFs + build full clickable PDF Outline Bookmark tree!
    print("=== PASS 4: Merging into single seamless IEEE PDF with full clickable bookmarks ===")
    writer = PdfWriter()
    for ppath in final_pdfs:
        writer.append(str(ppath), import_outline=False)

    total_pages = len(writer.pages)

    # Add Frontmatter & hierarchical Part I..IX -> Algorithm/Comparison bookmarks
    writer.add_outline_item("Frontmatter & Table of Contents", 0)
    cur_part_bm = None
    for entry in final_toc:
        pg_idx = max(0, min(total_pages - 1, entry["page"] - 1))
        if entry["level"] == 1:
            cur_part_bm = writer.add_outline_item(entry["title"], pg_idx)
        else:
            writer.add_outline_item(entry["title"], pg_idx, parent=cur_part_bm)

    writer.add_metadata({
        "/Title": "Mathematical Foundations, Formal Proofs, and Applications of 200 Clustering Algorithms (IEEE Monograph)",
        "/Author": "Gemini, AI made by Google",
        "/Subject": "200 Clustering Algorithms — Original Mathematical Proofs, Verified Classical Derivations, Applications, and Cross-Paradigm Comparisons",
        "/Creator": "Typst 0.15.1 + ClusterOverview Build Pipeline (Created by Gemini, AI made by Google)",
    })

    out_pdf = MONO / "ClusterOverview_200_Clustering_Algorithms_IEEE_Monograph.pdf"
    with open(out_pdf, "wb") as f:
        writer.write(f)

    size_mb = round(out_pdf.stat().st_size / (1024 * 1024), 2)
    print(f"\nSUCCESS!! Single Seamless IEEE Monograph compiled:")
    print(f"  File:  {out_pdf}")
    print(f"  Pages: {total_pages}")
    print(f"  Size:  {size_mb} MB")

    # Render 8 representative QA PNGs directly from the exact compiled volumes (low RAM, instant!)
    print("\n=== PASS 5: Rendering 8 representative QA PNG pages across the entire monograph ===")
    qa_targets = [
        (TMP_DIR / "vol0_frontmatter.typ", "1", "qa_01_title.png"),
        (TMP_DIR / "vol0_frontmatter.typ", "3", "qa_02_toc.png"),
        (TMP_DIR / "vol1_parts1_2.typ", "1", "qa_03_part1_opener.png"),
        (TMP_DIR / "vol2_parts3_4.typ", "105", "qa_04_cat4_algo.png"),
        (TMP_DIR / "vol4_parts7_8.typ", "120", "qa_05_cat8_algo.png"),
        (TMP_DIR / "vol5_part9_bib.typ", "1", "qa_06_part9_opener.png"),
        (TMP_DIR / "vol5_part9_bib.typ", "17", "qa_07_cmp2_equivalences.png"),
        (TMP_DIR / "vol5_part9_bib.typ", "172", "qa_08_cmp13_benchmarks.png"),
    ]
    for tpath, pg, out_name in qa_targets:
        out_png = QA_DIR / out_name
        r = subprocess.run(
            ["typst", "compile", "--root", str(MONO), "--format", "png", "--ppi", "110", "--pages", pg, str(tpath), str(out_png)],
            cwd=str(MONO),
            capture_output=True,
            text=True,
        )
        if r.returncode == 0 and out_png.exists():
            print(f"  Rendered {out_name} ({out_png.stat().st_size} bytes)")
        else:
            print(f"  Note on {out_name}: {r.stderr[:200]}")

    # Collect formal mathematical statistics across all 214 chapters
    algs = sorted(CH_DIR.glob("alg_*.typ"))
    cmps = sorted(CH_DIR.glob("cmp_*.typ"))
    stats = {
        "total_pdf_pages": total_pages,
        "pdf_size_mb": size_mb,
        "frontmatter_pages": fm_pages,
        "alg_chapters": len(algs),
        "cmp_chapters": len(cmps),
        "total_source_mb": round(sum(f.stat().st_size for f in algs + cmps) / (1024 * 1024), 2),
        "definitions": 0, "lemmas": 0, "theorems": 0, "propositions": 0, "corollaries": 0,
        "proofs": 0, "algorithms_pseudocode": 0, "tables": 0,
        "original_proof_badges": 0, "classical_cited_badges": 0,
        "bookmarks_count": len(final_toc) + 1,
    }
    for f in algs + cmps:
        t = f.read_text(encoding="utf-8")
        stats["definitions"] += len(re.findall(r"#definition\(", t))
        stats["lemmas"] += len(re.findall(r"#lemma\(", t))
        stats["theorems"] += len(re.findall(r"#theorem\(", t))
        stats["propositions"] += len(re.findall(r"#proposition\(", t))
        stats["corollaries"] += len(re.findall(r"#corollary\(", t))
        stats["proofs"] += len(re.findall(r"#proof(?:\[|\()", t))
        stats["algorithms_pseudocode"] += len(re.findall(r"#algbox\(", t))
        stats["tables"] += len(re.findall(r"table\(", t))
        stats["original_proof_badges"] += len(re.findall(r'origin:\s*"[^"]*(?:Original|Unified|Constructive)', t, re.I))
        stats["classical_cited_badges"] += len(re.findall(r'origin:\s*"[^"]*(?:Classical|Cited|Adapted)', t, re.I))
    stats["total_formal_statements"] = stats["lemmas"] + stats["theorems"] + stats["propositions"] + stats["corollaries"]

    (MONO / "data" / "monograph_stats.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print("\n=== MONOGRAPH STATS ===")
    print(json.dumps(stats, indent=2))

    # Clean up temporary volume PDFs to keep disk clean (keep .typ build script reproducible)
    for p in TMP_DIR.glob("*.pdf"):
        p.unlink()


if __name__ == "__main__":
    main()
