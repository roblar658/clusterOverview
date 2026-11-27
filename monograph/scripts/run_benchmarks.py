import csv
import json
import math
import os
import random
import subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(r"C:\dev\clusterOverview")
MONO = ROOT / "monograph"
DATA_DIR = MONO / "data"
FIG_DIR = MONO / "figures"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

CAT_NAMES = {
    1: "I. Centroid/Medoid",
    2: "II. Density/Mode",
    3: "III. Hierarchical",
    4: "IV. Distribution/Topic",
    5: "V. Graph/Spectral",
    6: "VI. Grid/Subspace",
    7: "VII. Fuzzy/Domain",
    8: "VIII. Neural/Stream/Evo",
}

CAT_COLORS = {
    1: "#1d4ed8",
    2: "#059669",
    3: "#7c3aed",
    4: "#d97706",
    5: "#dc2626",
    6: "#0891b2",
    7: "#be185d",
    8: "#475569",
}


def generate_synthetic_datasets():
    rng = random.Random(2026)

    # 1. Multi-shape non-convex dataset (2 moons + compact blob + elongated ellipse + 8 noise points)
    ms_path = DATA_DIR / "multishape_demo.csv"
    rows = []
    # Moon 1 (48 pts)
    for i in range(48):
        t = math.pi * i / 47.0
        x = math.cos(t) * 2.2 + rng.gauss(0, 0.12)
        y = math.sin(t) * 1.4 + rng.gauss(0, 0.12)
        rows.append((x, y, "moon_upper"))
    # Moon 2 (48 pts)
    for i in range(48):
        t = math.pi * i / 47.0
        x = 1.6 - math.cos(t) * 2.2 + rng.gauss(0, 0.12)
        y = 0.55 - math.sin(t) * 1.4 + rng.gauss(0, 0.12)
        rows.append((x, y, "moon_lower"))
    # Compact Gaussian blob (42 pts)
    for _ in range(42):
        x = -3.4 + rng.gauss(0, 0.35)
        y = -1.6 + rng.gauss(0, 0.35)
        rows.append((x, y, "compact_blob"))
    # Elongated rotated ellipse (42 pts)
    ang = math.radians(38)
    ca, sa = math.cos(ang), math.sin(ang)
    for _ in range(42):
        u = rng.gauss(0, 1.35)
        v = rng.gauss(0, 0.18)
        x = 4.6 + ca * u - sa * v
        y = -1.4 + sa * u + ca * v
        rows.append((x, y, "elongated_ellipse"))
    # Uniform outliers (8 pts)
    for _ in range(8):
        x = rng.uniform(-5.5, 6.8)
        y = rng.uniform(-3.6, 3.2)
        rows.append((x, y, "outlier_bg"))
    rng.shuffle(rows)
    with open(ms_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x1", "x2", "true_shape"])
        for x, y, lbl in rows:
            w.writerow([f"{x:.4f}", f"{y:.4f}", lbl])

    # 2. High-D anisotropic subspace + noise features + outliers (n=180, d=12, k=4)
    hd_path = DATA_DIR / "highdim_noisy_demo.csv"
    centers = [
        [2.8, 2.5, -2.0, 1.5],
        [-2.8, 2.2, 2.4, -1.5],
        [2.4, -2.6, 2.2, 2.0],
        [-2.5, -2.4, -2.5, -2.2],
    ]
    hd_rows = []
    for cid, c in enumerate(centers):
        for _ in range(42):
            sig = [c[j] + rng.gauss(0, 0.45 + 0.15 * j) for j in range(4)]
            noise = [rng.gauss(0, 1.8) for _ in range(8)]
            hd_rows.append((sig + noise, f"subspace_{cid+1}"))
    for _ in range(12):
        sig = [rng.uniform(-6.5, 6.5) for _ in range(4)]
        noise = [rng.uniform(-5.0, 5.0) for _ in range(8)]
        hd_rows.append((sig + noise, "outlier_hd"))
    rng.shuffle(hd_rows)
    with open(hd_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"s{j+1}" for j in range(4)] + [f"n{j+1}" for j in range(8)] + ["true_group"])
        for vals, lbl in hd_rows:
            w.writerow([f"{v:.4f}" for v in vals] + [lbl])

    return ms_path, hd_path


def run_rust_benchmarks(ms_path, hd_path):
    exe = ROOT / "target" / "debug" / "cluster_overview.exe"
    datasets = [
        ("iris", ROOT / "sample_data" / "iris_demo.csv", 3, "Iris (Convex/Overlapping, n=120, d=4)"),
        ("customers", ROOT / "sample_data" / "customers_mixed.csv", 4, "Customers Mixed (Num+Cat, n=150, d=6)"),
        ("multishape", ms_path, 4, "Multi-Shape (Non-Convex+Outliers, n=188, d=2)"),
        ("highdim", hd_path, 4, "High-D Subspace+Noise (n=180, d=12)"),
    ]
    all_runs = {}
    for key, path, k, label in datasets:
        out_csv = DATA_DIR / f"bench_{key}.csv"
        cmd = [
            str(exe),
            "--csv", str(path),
            "--metrics-csv", str(out_csv),
            "--skip-pdf",
            "--k", str(k),
            "--max-samples", "240",
        ]
        res = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
        print(f"[{key}] exit={res.returncode} :: {res.stdout.strip().splitlines()[-1] if res.stdout.strip() else ''}")
        if res.returncode != 0:
            raise RuntimeError(f"Benchmark {key} failed:\n{res.stderr}")
        rows = []
        with open(out_csv, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                rows.append({
                    "id": int(row["id"]),
                    "cat_num": int(row["cat_num"]),
                    "name": row["name"],
                    "n_clusters": int(row["n_clusters"]),
                    "n_noise": int(row["n_noise"]),
                    "noise_ratio": float(row["noise_ratio"]),
                    "silhouette": float(row["silhouette"]) if row["silhouette"] else None,
                    "davies_bouldin": float(row["davies_bouldin"]) if row["davies_bouldin"] else None,
                    "calinski_harabasz": float(row["calinski_harabasz"]) if row["calinski_harabasz"] else None,
                    "dunn_index": float(row["dunn_index"]) if row["dunn_index"] else None,
                    "inertia_sse": float(row["inertia_sse"]),
                    "balance_entropy": float(row["balance_entropy"]),
                    "ari": float(row["ari"]) if row["ari"] else None,
                    "duration_ms": float(row["duration_ms"]),
                })
        all_runs[key] = {"label": label, "k": k, "rows": rows}
    return all_runs


def summarize_and_plot(all_runs):
    ds_keys = ["iris", "customers", "multishape", "highdim"]
    summary = {"datasets": {}, "categories": {}, "top_overall": [], "per_algo": {}}

    # Per-algorithm aggregate across 4 datasets
    by_id = defaultdict(list)
    for key in ds_keys:
        rows = all_runs[key]["rows"]
        sils = [r["silhouette"] for r in rows if r["silhouette"] is not None]
        aris = [r["ari"] for r in rows if r["ari"] is not None]
        times = [r["duration_ms"] for r in rows]
        summary["datasets"][key] = {
            "label": all_runs[key]["label"],
            "mean_sil": sum(sils) / max(len(sils), 1),
            "max_sil": max(sils) if sils else 0.0,
            "mean_ari": sum(aris) / max(len(aris), 1) if aris else None,
            "max_ari": max(aris) if aris else None,
            "total_ms": sum(times),
            "median_ms": sorted(times)[len(times) // 2],
            "top5_sil": sorted(
                [r for r in rows if r["silhouette"] is not None],
                key=lambda x: x["silhouette"],
                reverse=True,
            )[:5],
            "top5_ari": sorted(
                [r for r in rows if r["ari"] is not None],
                key=lambda x: x["ari"],
                reverse=True,
            )[:5],
        }
        for r in rows:
            by_id[r["id"]].append((key, r))

    # Category summaries per dataset
    for cat in range(1, 9):
        cat_info = {"name": CAT_NAMES[cat], "datasets": {}}
        for key in ds_keys:
            crows = [r for r in all_runs[key]["rows"] if r["cat_num"] == cat]
            sils = [r["silhouette"] for r in crows if r["silhouette"] is not None]
            aris = [r["ari"] for r in crows if r["ari"] is not None]
            times = [r["duration_ms"] for r in crows]
            cat_info["datasets"][key] = {
                "count": len(crows),
                "mean_sil": round(sum(sils) / max(len(sils), 1), 4) if sils else 0.0,
                "max_sil": round(max(sils), 4) if sils else 0.0,
                "mean_ari": round(sum(aris) / max(len(aris), 1), 4) if aris else None,
                "max_ari": round(max(aris), 4) if aris else None,
                "mean_ms": round(sum(times) / max(len(times), 1), 3),
                "median_ms": round(sorted(times)[len(times) // 2], 3),
            }
        summary["categories"][str(cat)] = cat_info

    # Per-algorithm cross-dataset table
    algo_agg = []
    for aid in range(1, 201):
        entries = dict(by_id[aid])
        name = entries["iris"]["name"]
        cat = entries["iris"]["cat_num"]
        sils = [entries[k]["silhouette"] for k in ds_keys if entries[k]["silhouette"] is not None]
        aris = [entries[k]["ari"] for k in ds_keys if entries[k]["ari"] is not None]
        times = [entries[k]["duration_ms"] for k in ds_keys]
        rec = {
            "id": aid,
            "cat": cat,
            "name": name,
            "mean_sil": round(sum(sils) / max(len(sils), 1), 4) if sils else None,
            "mean_ari": round(sum(aris) / max(len(aris), 1), 4) if aris else None,
            "mean_ms": round(sum(times) / len(times), 3),
            "iris_sil": entries["iris"]["silhouette"],
            "iris_ari": entries["iris"]["ari"],
            "cust_sil": entries["customers"]["silhouette"],
            "ms_sil": entries["multishape"]["silhouette"],
            "ms_ari": entries["multishape"]["ari"],
            "hd_sil": entries["highdim"]["silhouette"],
            "hd_ari": entries["highdim"]["ari"],
        }
        algo_agg.append(rec)
        summary["per_algo"][str(aid)] = rec

    summary["top_overall"] = sorted(
        [a for a in algo_agg if a["mean_sil"] is not None],
        key=lambda x: (x["mean_sil"] + (x["mean_ari"] or 0.0)) / 2.0,
        reverse=True,
    )[:25]

    with open(DATA_DIR / "benchmark_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # --- SVG 1: Grouped Category Silhouette across 4 datasets ---
    svg_w, svg_h = 720, 320
    pad_l, pad_r, pad_t, pad_b = 56, 20, 36, 78
    pw = svg_w - pad_l - pad_r
    ph = svg_h - pad_t - pad_b
    ds_colors = {"iris": "#1d4ed8", "customers": "#059669", "multishape": "#d97706", "highdim": "#7c3aed"}
    ds_short = {"iris": "Iris", "customers": "Customers", "multishape": "Multi-Shape", "highdim": "High-D Noisy"}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="360" y="20" text-anchor="middle" font-family="Times New Roman, serif" font-size="13" font-weight="bold" fill="#0f172a">Mean Silhouette Score by Algorithm Category Across Four Benchmark Geometries (N = 200 Algorithms)</text>',
    ]
    # Y-axis grid 0.0 to 0.7
    ymax = 0.70
    for tick in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:
        y = pad_t + ph * (1.0 - tick / ymax)
        lines.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{svg_w-pad_r}" y2="{y:.1f}" stroke="#e2e8f0" stroke-width="0.8"/>')
        lines.append(f'<text x="{pad_l-6}" y="{y+3.5:.1f}" text-anchor="end" font-family="Times New Roman, serif" font-size="10" fill="#475569">{tick:.1f}</text>')
    lines.append(f'<line x1="{pad_l}" y1="{pad_t}" x2="{pad_l}" y2="{pad_t+ph}" stroke="#334155" stroke-width="1.1"/>')
    lines.append(f'<line x1="{pad_l}" y1="{pad_t+ph}" x2="{svg_w-pad_r}" y2="{pad_t+ph}" stroke="#334155" stroke-width="1.1"/>')

    gw = pw / 8.0
    bw = (gw - 16.0) / 4.0
    for cat in range(1, 9):
        gx = pad_l + (cat - 1) * gw + 8.0
        for j, k in enumerate(ds_keys):
            val = max(0.0, min(ymax, summary["categories"][str(cat)]["datasets"][k]["mean_sil"]))
            bh = ph * (val / ymax)
            x = gx + j * bw
            y = pad_t + ph - bh
            lines.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw-1.5:.1f}" height="{bh:.1f}" fill="{ds_colors[k]}" rx="1.5"/>')
        cx = gx + (gw - 16.0) / 2.0
        lbl = CAT_NAMES[cat].split(" ", 1)[0] + " " + CAT_NAMES[cat].split(" ", 1)[1].split("/")[0]
        lines.append(f'<text x="{cx:.1f}" y="{pad_t+ph+14}" text-anchor="middle" font-family="Times New Roman, serif" font-size="9.5" fill="#0f172a">{lbl}</text>')
    # Legend
    lx = 155
    for k in ds_keys:
        lines.append(f'<rect x="{lx}" y="{svg_h-26}" width="11" height="11" fill="{ds_colors[k]}" rx="1.5"/>')
        lines.append(f'<text x="{lx+15}" y="{svg_h-17}" font-family="Times New Roman, serif" font-size="10.5" fill="#1e293b">{ds_short[k]}</text>')
        lx += 110
    lines.append('</svg>')
    (FIG_DIR / "fig_cat_silhouette.svg").write_text("\n".join(lines), encoding="utf-8")

    # --- SVG 2: Log Runtime vs Mean Silhouette Scatter (200 Algorithms) ---
    svg_w, svg_h = 720, 340
    pad_l, pad_r, pad_t, pad_b = 58, 145, 36, 48
    pw = svg_w - pad_l - pad_r
    ph = svg_h - pad_t - pad_b
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="300" y="20" text-anchor="middle" font-family="Times New Roman, serif" font-size="13" font-weight="bold" fill="#0f172a">Quality–Speed Trade-off Across All 200 Clustering Algorithms (Mean Silhouette vs. Log10 Runtime)</text>',
    ]
    xmin, xmax = -1.5, 2.5  # log10(ms): 0.03 ms to 316 ms
    ymin, ymax = -0.1, 0.70
    for xt in [-1, 0, 1, 2]:
        x = pad_l + pw * ((xt - xmin) / (xmax - xmin))
        lines.append(f'<line x1="{x:.1f}" y1="{pad_t}" x2="{x:.1f}" y2="{pad_t+ph}" stroke="#e2e8f0" stroke-width="0.8"/>')
        lbl = f"{10**xt:g} ms"
        lines.append(f'<text x="{x:.1f}" y="{pad_t+ph+15}" text-anchor="middle" font-family="Times New Roman, serif" font-size="10" fill="#475569">{lbl}</text>')
    for yt in [0.0, 0.2, 0.4, 0.6]:
        y = pad_t + ph * (1.0 - (yt - ymin) / (ymax - ymin))
        lines.append(f'<line x1="{pad_l}" y1="{y:.1f}" x2="{pad_l+pw}" y2="{y:.1f}" stroke="#e2e8f0" stroke-width="0.8"/>')
        lines.append(f'<text x="{pad_l-6}" y="{y+3.5:.1f}" text-anchor="end" font-family="Times New Roman, serif" font-size="10" fill="#475569">{yt:.1f}</text>')
    lines.append(f'<rect x="{pad_l}" y="{pad_t}" width="{pw}" height="{ph}" fill="none" stroke="#334155" stroke-width="1.1"/>')
    lines.append(f'<text x="{pad_l+pw/2:.1f}" y="{svg_h-10}" text-anchor="middle" font-family="Times New Roman, serif" font-size="11" fill="#0f172a">Mean Wall-Clock Runtime per Dataset (Log Scale, ms)</text>')
    lines.append(f'<text x="16" y="{pad_t+ph/2:.1f}" text-anchor="middle" transform="rotate(-90 16 {pad_t+ph/2:.1f})" font-family="Times New Roman, serif" font-size="11" fill="#0f172a">Mean Silhouette Score</text>')
    for a in algo_agg:
        if a["mean_sil"] is None:
            continue
        lx = max(xmin, min(xmax, math.log10(max(a["mean_ms"], 0.03))))
        sy = max(ymin, min(ymax, a["mean_sil"]))
        px = pad_l + pw * ((lx - xmin) / (xmax - xmin))
        py = pad_t + ph * (1.0 - (sy - ymin) / (ymax - ymin))
        c = CAT_COLORS[a["cat"]]
        lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="3.4" fill="{c}" fill-opacity="0.78" stroke="#ffffff" stroke-width="0.5"/>')
    # Right legend
    ly = pad_t + 12
    for cat in range(1, 9):
        lines.append(f'<circle cx="{pad_l+pw+16}" cy="{ly}" r="4.2" fill="{CAT_COLORS[cat]}"/>')
        lines.append(f'<text x="{pad_l+pw+25}" y="{ly+3.5}" font-family="Times New Roman, serif" font-size="9.5" fill="#1e293b">{CAT_NAMES[cat]}</text>')
        ly += 22
    lines.append('</svg>')
    (FIG_DIR / "fig_ari_vs_runtime.svg").write_text("\n".join(lines), encoding="utf-8")

    # --- SVG 3: Category ARI Heatmap Across Ground-Truth Benchmarks ---
    svg_w, svg_h = 720, 290
    gt_keys = ["iris", "multishape", "highdim"]
    gt_labels = ["Iris (Convex)", "Multi-Shape (Non-Convex)", "High-D Subspace+Noise"]
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="360" y="20" text-anchor="middle" font-family="Times New Roman, serif" font-size="13" font-weight="bold" fill="#0f172a">Category Best &amp; Mean Adjusted Rand Index (ARI) vs. Ground Truth Across Diverse Geometries</text>',
    ]
    x0, y0, cw, ch = 178, 52, 172, 26
    for j, lbl in enumerate(gt_labels):
        lines.append(f'<text x="{x0 + j*cw + cw/2:.1f}" y="{y0-10}" text-anchor="middle" font-family="Times New Roman, serif" font-size="10.5" font-weight="bold" fill="#0f172a">{lbl}</text>')
    for cat in range(1, 9):
        y = y0 + (cat - 1) * ch
        lines.append(f'<text x="{x0-8}" y="{y+17}" text-anchor="end" font-family="Times New Roman, serif" font-size="10.5" fill="#0f172a">{CAT_NAMES[cat]}</text>')
        for j, k in enumerate(gt_keys):
            d = summary["categories"][str(cat)]["datasets"][k]
            ma = d["mean_ari"] or 0.0
            xa = d["max_ari"] or 0.0
            # Color intensity based on max ARI
            t = max(0.0, min(1.0, xa))
            r = int(241 - t * (241 - 29))
            g = int(245 - t * (245 - 78))
            b = int(249 - t * (249 - 216))
            tc = "#ffffff" if t > 0.58 else "#0f172a"
            x = x0 + j * cw
            lines.append(f'<rect x="{x+2}" y="{y+2}" width="{cw-4}" height="{ch-4}" fill="rgb({r},{g},{b})" rx="3" stroke="#cbd5e1" stroke-width="0.6"/>')
            lines.append(f'<text x="{x+cw/2:.1f}" y="{y+17}" text-anchor="middle" font-family="Times New Roman, serif" font-size="10" fill="{tc}">Mean {ma:.3f}  |  Max {xa:.3f}</text>')
    lines.append('</svg>')
    (FIG_DIR / "fig_geometry_heatmap.svg").write_text("\n".join(lines), encoding="utf-8")

    # --- SVG 4: Unified Mathematical Taxonomy Diagram ---
    svg_w, svg_h = 720, 310
    nodes = [
        (1, 95, 78, "Cat I: Centroid / Medoid (1-30)", "Bregman / Voronoi / L_p"),
        (4, 275, 78, "Cat IV: Distribution / Topic (81-102)", "EM / Variational / Exp-Family"),
        (7, 455, 78, "Cat VII: Fuzzy & Domain (151-175)", "Soft Membership / DTW / Ensembles"),
        (8, 630, 78, "Cat VIII: Neural & Stream (176-200)", "Deep Latent / Coresets / Meta"),
        (2, 95, 228, "Cat II: Density & Mode (31-55)", "KDE Level-Sets / Gradient Flow"),
        (3, 275, 228, "Cat III: Hierarchical (56-80)", "Ultrametrics / MST / Lance-Williams"),
        (5, 455, 228, "Cat V: Graph & Spectral (103-132)", "Laplacian Cuts / Modularity / MDL"),
        (6, 630, 228, "Cat VI: Grid & Subspace (133-150)", "Cell Quantization / Sparse L1"),
    ]
    edges = [
        (1, 4, "sigma->0"),
        (1, 7, "m->1"),
        (1, 5, "Kernel=Cut"),
        (2, 3, "MST Level-Set"),
        (2, 4, "KDE=GMM"),
        (3, 5, "Cut Tree"),
        (5, 6, "SSC Graph"),
        (4, 8, "VaDE/DEC"),
        (1, 8, "StreamKM++"),
        (7, 5, "Multi-View"),
    ]
    pos = {n[0]: (n[1], n[2]) for n in nodes}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_w}" height="{svg_h}" viewBox="0 0 {svg_w} {svg_h}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="360" y="20" text-anchor="middle" font-family="Times New Roman, serif" font-size="13" font-weight="bold" fill="#0f172a">Mathematical Unification Map Across All Eight Clustering Families (200 Algorithms)</text>',
    ]
    for u, v, lbl in edges:
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
        lines.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#94a3b8" stroke-width="1.2" stroke-dasharray="4,2"/>')
        lines.append(f'<rect x="{mx-31:.1f}" y="{my-7.5:.1f}" width="62" height="14" fill="#ffffff" stroke="#cbd5e1" stroke-width="0.6" rx="2"/>')
        lines.append(f'<text x="{mx:.1f}" y="{my+3:.1f}" text-anchor="middle" font-family="Times New Roman, serif" font-size="8.2" fill="#334155">{lbl}</text>')
    for cid, x, y, title, sub in nodes:
        c = CAT_COLORS[cid]
        lines.append(f'<rect x="{x-82}" y="{y-28}" width="164" height="56" fill="#f8fafc" stroke="{c}" stroke-width="1.8" rx="6"/>')
        lines.append(f'<rect x="{x-82}" y="{y-28}" width="164" height="7" fill="{c}" rx="3"/>')
        lines.append(f'<text x="{x}" y="{y-4}" text-anchor="middle" font-family="Times New Roman, serif" font-size="9.5" font-weight="bold" fill="#0f172a">{title}</text>')
        lines.append(f'<text x="{x}" y="{y+14}" text-anchor="middle" font-family="Times New Roman, serif" font-size="8.8" fill="#475569">{sub}</text>')
    lines.append('<text x="360" y="294" text-anchor="middle" font-family="Times New Roman, serif" font-size="9.5" font-style="italic" fill="#475569">Dashed edges denote formal mathematical equivalences, asymptotic limits, or exact reductions proved in Part IX (Chapter II).</text>')
    lines.append('</svg>')
    (FIG_DIR / "fig_taxonomy_map.svg").write_text("\n".join(lines), encoding="utf-8")

    print("Saved benchmark_summary.json and 4 SVG figures!")


if __name__ == "__main__":
    ms, hd = generate_synthetic_datasets()
    runs = run_rust_benchmarks(ms, hd)
    summarize_and_plot(runs)
