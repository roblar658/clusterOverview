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
per_algo_bench = bench.get("per_algo", {})

start_id = int(args.get("start_id", 1))
end_id = int(args.get("end_id", 200))
batch_size = int(args.get("batch_size", 4))
wave_size = int(args.get("wave_size", 5))
model_tier = args.get("model", "inherit")

BIB_KEYS_STR = ", ".join(bib_keys)

CAT_TITLES = {
    1: "Category I: Centroid, Medoid & Partitioning (Algs 1-30)",
    2: "Category II: Density & Mode-Seeking (Algs 31-55)",
    3: "Category III: Hierarchical & Tree-Based (Algs 56-80)",
    4: "Category IV: Distribution, Statistical & Topic Models (Algs 81-102)",
    5: "Category V: Graph, Spectral & Network (Algs 103-132)",
    6: "Category VI: Grid, Subspace & Projection (Algs 133-150)",
    7: "Category VII: Fuzzy, Biclustering, Time-Series & Domain (Algs 151-175)",
    8: "Category VIII: Neural, Deep, Stream & Evolutionary (Algs 176-200)",
}

quota_hit = False


def verify_ids_on_disk(ids: list) -> dict:
    """Batch-verify multiple algorithm IDs in a single Python invocation."""
    existing = []
    out = {}
    for aid in ids:
        f = CH_DIR / f"alg_{aid:03d}.typ"
        if not f.exists() or f.stat().st_size < 5000:
            out[aid] = (False, 0, "missing_or_small")
        else:
            existing.append(aid)
    if not existing:
        return out
    rel_files = [f"chapters/alg_{aid:03d}.typ" for aid in existing]
    res = subprocess.run(
        ["python", str(MONO / "scripts" / "verify_chapter.py"), "--range", "4", "5", *rel_files],
        cwd=str(MONO),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    try:
        parsed = json.loads(res.stdout)
        for aid, rec in zip(existing, parsed):
            ok = bool(rec.get("ok") and rec.get("length_ok") and not rec.get("warnings"))
            pgs = int(rec.get("content_pages", 0))
            err = ""
            if not rec.get("ok"):
                err = (rec.get("stderr") or rec.get("error") or "compile_error")[:400]
            elif not rec.get("length_ok"):
                err = f"pages={pgs} not in [4,5]"
            elif rec.get("warnings"):
                err = f"warnings={rec['warnings'][:2]}"
            out[aid] = (ok, pgs, err)
    except Exception as e:
        for aid in existing:
            out[aid] = (False, 0, f"parse_err: {e}")
    return out


def format_algo_brief(item: dict, existing_note: str = "") -> str:
    aid = item["id"]
    b = per_algo_bench.get(str(aid), {})
    def fmt(v):
        return f"{v:.3f}" if isinstance(v, (int, float)) else "N/A"
    extra = f"\n  NOTE (existing file needs fix): {existing_note}" if existing_note else ""
    return (
        f"- Algorithm #{aid}: \"{item['name']}\" | Category: \"{item['category']}\"\n"
        f"  Summary: {item['desc_no']}\n"
        f"  Target file: C:\\dev\\clusterOverview\\monograph\\chapters\\alg_{aid:03d}.typ\n"
        f"  ClusterOverview 4-dataset benchmark: "
        f"Iris Sil={fmt(b.get('iris_sil'))}, ARI={fmt(b.get('iris_ari'))}; "
        f"Customers Sil={fmt(b.get('cust_sil'))}; "
        f"Multi-Shape Sil={fmt(b.get('ms_sil'))}, ARI={fmt(b.get('ms_ari'))}; "
        f"High-D Sil={fmt(b.get('hd_sil'))}, ARI={fmt(b.get('hd_ari'))}; "
        f"Mean runtime={fmt(b.get('mean_ms'))} ms.{extra}"
    )


def build_prompt(batch: list, notes_map: dict) -> str:
    briefs = "\n".join(format_algo_brief(it, notes_map.get(it["id"], "")) for it in batch)
    ids_str = ", ".join(str(it["id"]) for it in batch)
    rel_list = " ".join(f"chapters/alg_{it['id']:03d}.typ" for it in batch)
    return f"""You are a senior applied mathematician and IEEE Transactions author contributing dedicated 4-page technical chapters to the monograph "Mathematical Foundations, Formal Proofs, and Applications of 200 Clustering Algorithms" (Created by Gemini, AI made by Google).

TASK: Author and verify the following {len(batch)} individual IEEE two-column Typst 0.15 chapters (Algorithms {ids_str}):
{briefs}

DO NOT read other chapter files — everything you need for exact template compliance is right here! Work in maximum efficiency: write ALL {len(batch)} files in parallel in Turn 1 (or edit only broken files if a NOTE says an existing file needs a small fix), then verify ALL {len(batch)} in ONE command in Turn 2!

EXACT TEMPLATE HEADER & STRUCTURE (~225–265 lines / ~24–27 KB per file -> lands cleanly on EXACTLY 4 IEEE pages!):
```typst
#import "../template/ieee.typ": *

#alg_entry(
  id: <ID>,
  name: "<English IEEE Title>",
  category: "<Category Label>",
  year: "<Year>",
  paradigm: "<Short Mathematical Paradigm>",
  complexity: $<Exact Time Complexity>$,
  provenance: "Original Proofs + Cited Classical Bounds",
)
```
Required 5 IEEE sections per file:
- `#sec("I", "Problem Formulation and Mathematical Model")`: Formal objective/probabilistic/graph/topological formulation, `#definition(title: "...")[...]`, and `#algbox(title: "...", inputs: [...], outputs: [...])[...]` with numbered mathematical steps.
- `#sec("II", "Formal Mathematical Proofs and Theoretical Guarantees")`: Subsections `#subsec("A", ...)` .. `#subsec("D", ...)` (or `"E"`) containing **4 to 5 formal numbered statements** (`#lemma`, `#theorem`, `#proposition`, `#corollary`) EACH followed by a rigorous, complete, step-by-step `#proof[...]`.
  * At least 2–3 results MUST be genuine original / constructive derivations (`origin: "Original Proof"`, `"Original Derivation"`, `"Original Counterexample"`, or `"Original Unified Proof"`): Lyapunov/Bregman/MM/variational monotonicity, constructive finite-sample counterexamples, geometric/scale/permutation invariance & breakdown-point proofs, algebraic reductions, or tight complexity bounds.
  * Reproduced/adapted landmark theorems MUST use `origin: "Classical — Cited"` and cite real `@key` references from the verified list below.
- `#sec("III", "Complexity, Geometry and Operational Profile")`: Exact time/space/memory/cache analysis + `#figure(table(columns: (1.1fr, 1.9fr), align: (left, left), [Property], [Characterization], ...), caption: [...])`.
- `#sec("IV", "Engineering and Scientific Applications")`:
  * `#subsec("A", "Primary Application Domains and Mathematical Fit")`: 5 detailed real-world scientific/engineering applications explaining mathematically *why* this algorithm's inductive bias matches the domain.
  * `#subsec("B", "Failure Modes and Domain Suitability Matrix")`: Explicit failure regimes + `#figure(table(columns: (1.0fr, 1.2fr, 0.8fr), align: (left, left, left), [Domain], [Mathematical Fit], [Caution / Alternative], ...), caption: [...])`.
- `#sec("V", "Structural Relations and Empirical Behavior")`: Formal `#proposition` or tight mathematical link to adjacent algorithms in the 200-catalog + empirical discussion referencing the real ClusterOverview 4-dataset benchmark metrics given above!

CRITICAL TYPST 0.15 SYNTAX & ANTI-OVERFLOW RULES (ZERO ERRORS ON FIRST COMPILE!):
1. Bold is SINGLE `*bold*` (NEVER `**bold**`!). Italic is `_italic_`.
2. Inner products: `chevron.l u, v chevron.r` (NEVER `angle.l`/`angle.r`!).
3. Space before parens after subscript: `D_phi (x || y)`, `hat(f)_h (x)`, `V_j (bold(mu))`.
4. Literal hash: write `Alg. 12` or `\\#12` (NEVER bare `#12`!).
5. Quote ALL multi-letter math operators: `"Tr"`, `"Var"`, `"Cov"`, `"diag"`, `"KL"`, `"cut"`, `"vol"`, `"deg"`, `"rank"`, `"span"`, `"dist"`, `"med"`, `"sign"`, `"poly"`, `"supp"`, `"ELBO"`, `"BIC"`, `"MDL"`, `"ARI"`, `"SSE"`, `"SBD"`, `"DTW"`, `"Ncut"`, `"Rcut"`, `"kNN"`, `"RNN"`, `"SNN"`, `"CF"`, `"Gini"`, `"LOF"`, `"Pr"`. Use `arg min` and `arg max` (with space).
6. Set, order & calculus symbols in Typst 0.15: intersection is `inter` (NEVER `sect`!), union is `union`, big product is `product` (NEVER `prod`!), PSD order is `succ.eq` / `prec.eq` (NEVER `succeq`/`preceq`!), equivalence/distributed-as is `tilde.op` or `~` (NEVER `sim`!), Minkowski/direct sum is `plus.o` (NEVER `plus.circle`!), differential/gradient are `dif` / `nabla` / `partial` (NEVER `diff`!). NEVER use LaTeX backslash commands like `\\text{...}`, `\\epsilon`, `\\lambda`, `\\sum`, `\\frac` inside Typst math — use `"..."`, `epsilon`, `lambda`, `sum`, `a / b`!
7. Escape literal asterisks in prose: write `DBSCAN\\*` or `A\\*`, and ALWAYS put spaces around `/` between bold/italic words (`*a* / *b*`, NEVER `*a*/*b*` because `*/` triggers a Typst block-comment syntax error!).
8. Keep EVERY display equation <= 65 chars wide (split multi-term formulas across lines with `&` and `\\`) so NOTHING overflows the 3.5-inch IEEE column!
9. NEVER use raw `@` except for citations, and ONLY cite verified BibTeX keys from this exact list:
[{BIB_KEYS_STR}]

VERIFICATION (ONE SINGLE COMMAND):
Run:
`python C:\\dev\\clusterOverview\\monograph\\scripts\\verify_chapter.py --range 4 5 {rel_list}`
All {len(batch)} entries MUST report `"ok": true`, `"length_ok": true` (`content_pages` 4 or 5, ideally 4!), and `"warnings": []`. Fix any issue and re-verify before calling `finish`.
"""


BATCH_SCHEMA = {
    "type": "object",
    "properties": {
        "verified_ids": {"type": "array", "items": {"type": "integer"}},
        "page_counts": {"type": "object"},
        "notes": {"type": "string"},
    },
    "required": ["verified_ids", "page_counts"],
}


async def run_batch_safe(batch: list, attempt: int = 1):
    global quota_hit
    if quota_hit:
        return {"verified_ids": [], "page_counts": {}, "failed_ids": [it["id"] for it in batch]}

    status_map = verify_ids_on_disk([it["id"] for it in batch])
    todo = []
    already = {}
    notes_map = {}
    for it in batch:
        ok, pgs, err = status_map[it["id"]]
        if ok:
            already[str(it["id"])] = pgs
        else:
            todo.append(it)
            if err and err != "missing_or_small":
                notes_map[it["id"]] = err

    if not todo:
        log(f"SKIP batch {[it['id'] for it in batch]} — already verified: {already}")
        return {"verified_ids": [it["id"] for it in batch], "page_counts": already, "failed_ids": []}

    ids_label = "-".join(f"{it['id']:03d}" for it in todo)
    role = f"Math Author Algs {ids_label}" + (f" (retry {attempt})" if attempt > 1 else "")
    try:
        await agent(
            build_prompt(todo, notes_map),
            schema=BATCH_SCHEMA,
            role=role,
            type_name="self",
            workspace="inherit",
            model=model_tier,
        )
        post = verify_ids_on_disk([it["id"] for it in todo])
        still_bad = []
        counts = dict(already)
        for it in todo:
            ok, pgs, err = post[it["id"]]
            if ok:
                counts[str(it["id"])] = pgs
            else:
                still_bad.append(it)
                log(f"  POST-CHECK FAILED alg_{it['id']:03d}: {err}")
        if still_bad and attempt < 2 and not quota_hit:
            log(f"RETRYING {[it['id'] for it in still_bad]} (attempt {attempt+1})")
            sub = await run_batch_safe(still_bad, attempt=attempt + 1)
            counts.update(sub["page_counts"])
            return {
                "verified_ids": [int(k) for k in counts.keys()],
                "page_counts": counts,
                "failed_ids": sub["failed_ids"],
            }
        if not still_bad:
            log(f"VERIFIED batch {[it['id'] for it in batch]} -> pages {counts}")
        return {
            "verified_ids": [int(k) for k in counts.keys()],
            "page_counts": counts,
            "failed_ids": [it["id"] for it in still_bad],
        }
    except Exception as exc:
        msg = str(exc)
        if "Individual quota reached" in msg:
            quota_hit = True
            log(f"INDIVIDUAL QUOTA HIT in batch {[it['id'] for it in todo]} — halting cleanly!")
        elif ("429" in msg or "ResourceExhausted" in msg) and attempt < 3:
            wait_s = 20 * attempt
            log(f"TRANSIENT 429 in batch {[it['id'] for it in todo]} (attempt {attempt}) — sleeping {wait_s}s and retrying...")
            await asyncio.sleep(wait_s)
            return await run_batch_safe(todo, attempt=attempt + 1)
        elif "429" in msg or "ResourceExhausted" in msg:
            quota_hit = True
            log(f"PERSISTENT 429 in batch {[it['id'] for it in todo]} after retries — halting cleanly!")
        else:
            log(f"ERROR in batch {[it['id'] for it in todo]} (attempt {attempt}): {msg[:300]}")
            if attempt < 2 and not quota_hit:
                return await run_batch_safe(todo, attempt=attempt + 1)
        # Re-check disk in case some files WERE written before error
        post = verify_ids_on_disk([it["id"] for it in todo])
        counts = dict(already)
        failed = []
        for it in todo:
            ok, pgs, _ = post[it["id"]]
            if ok:
                counts[str(it["id"])] = pgs
            else:
                failed.append(it["id"])
        return {
            "verified_ids": [int(k) for k in counts.keys()],
            "page_counts": counts,
            "failed_ids": failed,
        }


selected = [it for it in catalog if start_id <= it["id"] <= end_id]
by_cat = {}
for it in selected:
    by_cat.setdefault(it["cat_num"], []).append(it)

all_verified = {}
all_failed = []

for cat_num in sorted(by_cat.keys()):
    if quota_hit:
        all_failed.extend(it["id"] for it in by_cat[cat_num])
        continue
    phase(CAT_TITLES[cat_num])
    items = by_cat[cat_num]
    # Pre-check entire category in one fast batch call!
    pre = verify_ids_on_disk([it["id"] for it in items])
    todo_items = []
    for it in items:
        ok, pgs, _ = pre[it["id"]]
        if ok:
            all_verified[str(it["id"])] = pgs
        else:
            todo_items.append(it)
    if not todo_items:
        log(f"Category {cat_num}: ALL {len(items)} algorithms already verified on disk!")
        continue
    batches = [todo_items[i : i + batch_size] for i in range(0, len(todo_items), batch_size)]
    log(f"Category {cat_num}: {len(items)-len(todo_items)} already verified, {len(todo_items)} remaining in {len(batches)} batches")
    for w_idx in range(0, len(batches), wave_size):
        if quota_hit:
            for rem_b in batches[w_idx:]:
                all_failed.extend(it["id"] for it in rem_b)
            break
        wave = batches[w_idx : w_idx + wave_size]
        wave_results = await parallel(wave, lambda b: run_batch_safe(b))
        for r in wave_results:
            all_verified.update(r["page_counts"])
            all_failed.extend(r["failed_ids"])
    log(f"Category {cat_num} done: cumulative verified={len(all_verified)}, failed={len(all_failed)}, quota_hit={quota_hit}")

total_pages = sum(int(v) for v in all_verified.values())
log(f"WORKFLOW DONE: verified={len(all_verified)}/{len(selected)}, failed={len(all_failed)}, quota_hit={quota_hit}, total_pages={total_pages}")
print(json.dumps({
    "verified_count": len(all_verified),
    "failed_count": len(all_failed),
    "failed_ids": all_failed,
    "quota_hit": quota_hit,
    "total_pages": total_pages,
}, indent=2))
