import re
import subprocess
from pathlib import Path

MONO = Path(r"C:\dev\clusterOverview\monograph")
CH = MONO / "chapters"

# Scan alg_001..200 for any known stray LaTeX/Typst-v0.12 tokens in math
changed = []
for aid in range(1, 201):
    p = CH / f"alg_{aid:03d}.typ"
    if not p.exists():
        continue
    txt = p.read_text(encoding="utf-8")
    orig = txt
    # Only touch well-known invalid Typst 0.15 math identifiers
    txt = re.sub(r"\bnotin\b", "in.not", txt)
    txt = re.sub(r"\bcdot\b", "dot.c", txt)
    txt = re.sub(r"\bsucceq\b", "succ.eq", txt)
    txt = re.sub(r"\bpreceq\b", "prec.eq", txt)
    txt = re.sub(r"\\text\{([^}]*)\}", r'"\1"', txt)
    for gr in ["psi", "phi", "theta", "epsilon", "lambda", "sigma", "omega", "ell", "alpha", "beta", "gamma", "delta", "mu", "rho", "tau", "pi", "kappa", "eta", "zeta", "nu", "xi", "chi", "Gamma", "Delta", "Lambda", "Sigma", "Omega", "Phi", "Psi", "Theta"]:
        txt = re.sub(rf"\\{gr}(?![A-Za-z])", gr, txt)
    txt = re.sub(r"\\infty(?![A-Za-z])", "infinity", txt)
    if txt != orig:
        p.write_text(txt, encoding="utf-8")
        changed.append(aid)

print(f"Auto-cleaned {len(changed)} files: {changed}")
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
    print("STDERR:", res.stderr[:3000])
    raise SystemExit(res.returncode)
