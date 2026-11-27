import json
import re
import subprocess
from pathlib import Path

MONO = Path(__file__).resolve().parent.parent
CH_DIR = MONO / "chapters"
QA_DIR = MONO / "_qa_preview"
QA_DIR.mkdir(parents=True, exist_ok=True)

# 1. Build main.typ & compile full PDF
res = subprocess.run(
    ["python", str(MONO / "scripts" / "build_monograph.py"), "--compile"],
    cwd=str(MONO),
    capture_output=True,
    text=True,
    encoding="utf-8",
    errors="replace",
)
print(res.stdout)
if res.returncode != 0:
    print("STDERR:\n", res.stderr[:4000])
    raise SystemExit(res.returncode)

# 2. Collect formal mathematical statistics across all 214 chapters
algs = sorted(CH_DIR.glob("alg_*.typ"))
cmps = sorted(CH_DIR.glob("cmp_*.typ"))
all_files = algs + cmps
total_bytes = sum(f.stat().st_size for f in all_files)

stats = {
    "alg_chapters": len(algs),
    "cmp_chapters": len(cmps),
    "total_source_mb": round(total_bytes / (1024 * 1024), 2),
    "definitions": 0,
    "lemmas": 0,
    "theorems": 0,
    "propositions": 0,
    "corollaries": 0,
    "proofs": 0,
    "algorithms_pseudocode": 0,
    "tables": 0,
    "original_proof_badges": 0,
    "classical_cited_badges": 0,
}

for f in all_files:
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

stats["total_formal_statements"] = (
    stats["lemmas"] + stats["theorems"] + stats["propositions"] + stats["corollaries"]
)
print("MONOGRAPH STATS:")
print(json.dumps(stats, indent=2))

# 3. Query Part I..IX starting page numbers via typst query
q = subprocess.run(
    ["typst", "query", str(MONO / "main.typ"), "heading.where(level: 1)"],
    cwd=str(MONO),
    capture_output=True,
    text=True,
    encoding="utf-8",
    errors="replace",
)
if q.returncode == 0 and q.stdout.strip():
    h1s = json.loads(q.stdout)
    print(f"Found {len(h1s)} top-level Part/Section headings.")

# 4. Render 8 representative QA pages across the ~1030-page PDF using typst --pages
# We can pass --pages 1,10,430,835,860,990,1010,1030 (or exact ranges)
pdf_path = MONO / "ClusterOverview_200_Clustering_Algorithms_IEEE_Monograph.pdf"
pdf_bytes = pdf_path.read_bytes()
matches = re.findall(rb"/Type\s*/Pages\b[^>]*?/Count\s+(\d+)", pdf_bytes, flags=re.DOTALL)
total_pages = max(int(m) for m in matches) if matches else 0

sample_pages = sorted(set([
    1,                              # Title page + Abstract (Created by Gemini, AI made by Google)
    3,                              # Table of Contents
    max(1, total_pages // 4),       # Mid Part II/III Algorithm chapter
    max(1, total_pages // 2),       # Mid Part V Graph/Spectral chapter
    max(1, total_pages - 195),      # Part IX Comparison Ch I (Taxonomy SVG)
    max(1, total_pages - 175),      # Part IX Comparison Ch II (Equivalence proofs)
    max(1, total_pages - 35),       # Part IX Comparison Ch XIII (Empirical SVG charts)
    max(1, total_pages - 3),        # Complete IEEE Bibliography
]))
pages_arg = ",".join(str(p) for p in sample_pages)
print(f"Rendering QA sample pages: {pages_arg} (out of {total_pages} total pages)")
r_png = subprocess.run(
    [
        "typst", "compile", "--format", "png", "--ppi", "110",
        "--pages", pages_arg,
        str(MONO / "main.typ"),
        str(QA_DIR / "qa-{0p}.png"),
    ],
    cwd=str(MONO),
    capture_output=True,
    text=True,
    encoding="utf-8",
    errors="replace",
)
if r_png.returncode != 0:
    print("PNG render warning:", r_png.stderr[:1000])
else:
    for p in sorted(QA_DIR.glob("qa-*.png")):
        print(f"  QA PNG: {p.name} ({p.stat().st_size} bytes)")
