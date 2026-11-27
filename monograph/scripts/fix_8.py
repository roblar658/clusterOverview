import json
import re
import subprocess
from pathlib import Path

MONO = Path(r"C:\dev\clusterOverview\monograph")
CH = MONO / "chapters"

bad_ids = [36, 37, 40, 44, 45, 46, 47, 48]

# Common deterministic Typst 0.15 fixes across these 8 files
for aid in bad_ids:
    p = CH / f"alg_{aid:03d}.typ"
    txt = p.read_text(encoding="utf-8")
    # 1. Fix LaTeX \text{...} -> "..."
    txt = re.sub(r"\\text\{([^}]*)\}", r'"\1"', txt)
    # 2. Fix LaTeX \epsilon, \lambda, \sigma, \delta, \alpha, \beta, \gamma, \mu, \rho, \tau, \phi, \psi, \omega, \theta, \pi
    for gr in ["epsilon", "lambda", "sigma", "delta", "alpha", "beta", "gamma", "mu", "rho", "tau", "phi", "psi", "omega", "theta", "pi", "infty"]:
        repl = "infinity" if gr == "infty" else gr
        txt = re.sub(rf"\\{gr}\b", repl, txt)
    # 3. Fix Typst `prod_` / `prod ` -> `product_` / `product ` inside math (simple safe word boundary)
    txt = re.sub(r"\bprod(?=[_( ])", "product", txt)
    # 4. Fix `diff` in math (`diff ` -> `dif ` or `partial`) -> let's inspect first if still broken
    # 5. Fix unescaped `*` inside `$...$`? No, inside prose `($...$)` or `foo*` -> inspect unclosed delimiter lines first!
    p.write_text(txt, encoding="utf-8")

# Re-verify these 8 and print full diagnostics for any still failing
targets = [f"chapters/alg_{aid:03d}.typ" for aid in bad_ids]
res = subprocess.run(
    ["python", str(MONO / "scripts" / "verify_chapter.py"), "--range", "4", "5", *targets],
    cwd=str(MONO),
    capture_output=True,
    text=True,
    encoding="utf-8",
    errors="replace",
)
for r in json.loads(res.stdout):
    if r.get("ok") and r.get("length_ok") and not r.get("warnings"):
        print(f"OK: {r['file']} ({r['content_pages']} pages)")
    else:
        print(f"STILL BAD: {r['file']} pages={r.get('content_pages')} warn={r.get('warnings')}")
        if r.get("stderr"):
            print(r["stderr"][:1200])
