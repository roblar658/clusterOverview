import json
import subprocess
import sys
import tempfile
from pathlib import Path

MONO = Path(r"C:\dev\clusterOverview\monograph")


def verify_one(rel_path: str, min_pages: int = 4, max_pages: int = 5) -> dict:
    target = Path(rel_path)
    if not target.is_absolute():
        target = MONO / rel_path
    if not target.exists():
        return {"ok": False, "file": rel_path, "error": f"File not found: {target}"}

    rel_from_mono = target.relative_to(MONO).as_posix()
    wrapper_src = f'''#import "template/ieee.typ": *
#show: ieee_doc.with(title: "Verification Wrapper")
#include "{rel_from_mono}"
#context [#metadata(counter(page).get().first()) <chapter_end>]
#pagebreak()
#bibliography("references.bib", title: "References", style: "ieee")
'''
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        wrap_typ = MONO / f"_tmp_verify_{target.stem}.typ"
        out_pdf = td_path / "out.pdf"
        try:
            wrap_typ.write_text(wrapper_src, encoding="utf-8")
            res = subprocess.run(
                ["typst", "compile", str(wrap_typ), str(out_pdf)],
                cwd=str(MONO),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            if res.returncode != 0:
                return {
                    "ok": False,
                    "length_ok": False,
                    "file": rel_from_mono,
                    "returncode": res.returncode,
                    "stderr": res.stderr[:2500],
                }

            q = subprocess.run(
                ["typst", "query", str(wrap_typ), "<chapter_end>", "--field", "value", "--one"],
                cwd=str(MONO),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            pages = 0
            if q.returncode == 0 and q.stdout.strip():
                try:
                    pages = int(json.loads(q.stdout))
                except Exception:
                    pages = 0

            warnings = [ln for ln in res.stderr.splitlines() if ln.strip()]
            ok_len = (min_pages <= pages <= max_pages)
            return {
                "ok": True,
                "length_ok": ok_len,
                "file": rel_from_mono,
                "content_pages": pages,
                "target_range": [min_pages, max_pages],
                "warnings": warnings[:10],
            }
        finally:
            if wrap_typ.exists():
                wrap_typ.unlink()


if __name__ == "__main__":
    # Supports:
    #   python verify_chapter.py chapters/alg_001.typ 4 5
    #   python verify_chapter.py --range 4 5 chapters/alg_031.typ chapters/alg_032.typ ...
    argv = sys.argv[1:]
    if not argv:
        print("Usage: python verify_chapter.py [--range MIN MAX] <file1> [file2 ...]")
        sys.exit(2)

    if argv[0] == "--range":
        lo, hi = int(argv[1]), int(argv[2])
        files = argv[3:]
        results = [verify_one(f, lo, hi) for f in files]
        print(json.dumps(results, indent=2))
        all_ok = all(r.get("ok") and r.get("length_ok") and not r.get("warnings") for r in results)
        sys.exit(0 if all_ok else 1)
    else:
        f = argv[0]
        lo = int(argv[1]) if len(argv) > 1 else 4
        hi = int(argv[2]) if len(argv) > 2 else 5
        r = verify_one(f, lo, hi)
        print(json.dumps(r, indent=2))
        sys.exit(0 if r.get("ok") else 1)
