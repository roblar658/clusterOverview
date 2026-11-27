import json
import re
from pathlib import Path

src = Path(r"C:\dev\clusterOverview\src\catalog.rs").read_text(encoding="utf-8")
pat = re.compile(
    r'AlgorithmSpec\s*\{\s*id:\s*(\d+),\s*category:\s*Category::(\w+),\s*name:\s*"([^"]+)",\s*description_no:\s*"([^"]+)",?\s*\}',
    re.S,
)
cat_map = {
    "Centroid": (1, "Cat. I: Centroid & Partitioning"),
    "Density": (2, "Cat. II: Density & Mode-Seeking"),
    "Hierarchical": (3, "Cat. III: Hierarchical & Tree-Based"),
    "Distribution": (4, "Cat. IV: Distribution, Statistical & Topic"),
    "GraphSpectral": (5, "Cat. V: Graph, Spectral & Network"),
    "GridSubspace": (6, "Cat. VI: Grid, Subspace & Projection"),
    "FuzzyDomain": (7, "Cat. VII: Fuzzy, Biclustering & Domain"),
    "NeuralStream": (8, "Cat. VIII: Neural, Deep, Stream & Evolutionary"),
}
items = []
for m in pat.finditer(src):
    aid, cat_enum, name, desc = int(m.group(1)), m.group(2), m.group(3), m.group(4)
    cnum, clabel = cat_map[cat_enum]
    items.append({
        "id": aid,
        "cat_num": cnum,
        "cat_enum": cat_enum,
        "category": clabel,
        "name": name,
        "desc_no": desc,
    })
assert len(items) == 200, f"Expected 200, got {len(items)}"
Path(r"C:\dev\clusterOverview\monograph\data\catalog_200.json").write_text(
    json.dumps(items, indent=2, ensure_ascii=False), encoding="utf-8"
)

bib = Path(r"C:\dev\clusterOverview\monograph\references.bib").read_text(encoding="utf-8")
keys = re.findall(r"@\w+\{([^,\s]+),", bib)
Path(r"C:\dev\clusterOverview\monograph\data\bib_keys.json").write_text(
    json.dumps(sorted(set(keys)), indent=2), encoding="utf-8"
)
print(f"Exported {len(items)} algorithms and {len(set(keys))} verified BibTeX keys.")
