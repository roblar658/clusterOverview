import json
import re
import subprocess
from pathlib import Path

MONO = Path(r"C:\dev\clusterOverview\monograph")
CH = MONO / "chapters"

# Fix known Typst 0.15 syntax issues across existing chapters 31..48
for aid in range(31, 49):
    p = CH / f"alg_{aid:03d}.typ"
    if not p.exists():
        continue
    txt = p.read_text(encoding="utf-8")
    orig = txt
    # Replace deprecated Typst `sect` -> `inter` and `sect.big` -> `inter.big` inside math
    txt = re.sub(r"\bsect\.big\b", "inter.big", txt)
    txt = re.sub(r"\bsect\b", "inter", txt)
    # Fix `*word*/*word*` where `*/` is parsed as end-of-block-comment
    txt = txt.replace("*/*", "* / *")
    # Fix unescaped `DBSCAN*` in prose
    txt = re.sub(r"DBSCAN\*(?!\w)", r"DBSCAN\*", txt)
    if txt != orig:
        p.write_text(txt, encoding="utf-8")

targets = [f"chapters/alg_{aid:03d}.typ" for aid in range(1, 49)]
res = subprocess.run(
    ["python", str(MONO / "scripts" / "verify_chapter.py"), "--range", "4", "5", *targets],
    cwd=str(MONO),
    capture_output=True,
    text=True,
    encoding="utf-8",
    errors="replace",
)
data = json.loads(res.stdout)
ok_list = [r for r in data if r.get("ok") and r.get("length_ok") and not r.get("warnings")]
bad_list = [r for r in data if not (r.get("ok") and r.get("length_ok") and not r.get("warnings"))]
total_pages = sum(r.get("content_pages", 0) for r in ok_list)
print(f"VERIFIED: {len(ok_list)}/48 chapters ({total_pages} pages). NEED FIX: {len(bad_list)}")
for b in bad_list:
    print("BAD:", b["file"], "ok=", b.get("ok"), "pages=", b.get("content_pages"), "warn=", b.get("warnings"), "err=", (b.get("stderr") or "")[:300])
