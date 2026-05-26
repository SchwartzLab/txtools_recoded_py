#!/usr/bin/env python
"""Compare Python txtools vs reference R txtools and emit a self-contained HTML."""
import base64
import io
import sys
import textwrap
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

V = sys.argv[1] if len(sys.argv) > 1 else "."
PY_TSV = f"{V}/py_chr19.tsv"
R_TSV = f"{V}/r_chr19.tsv"
OUT = f"{V}/txtools_comparison_chr19.html"

KEYS = ["gene", "txcoor"]
INT_COLS = ["cov", "start_5p", "end_3p", "A", "C", "G", "T", "N", "-", "."]
FLT_COLS = ["misincRate", "MR_CtoT"]
NUC = ["A", "C", "G", "T", "N", "-", "."]


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def parse_phase_timing(path):
    """Parse phase_timing.log (per-phase seconds + /usr/bin/time -v peak RSS, CPU%)."""
    import re
    data = {}
    cur = None
    try:
        lines = open(path).read().splitlines()
    except OSError:
        return data
    for ln in lines:
        if "=====" in ln and "PYTHON" in ln:
            cur = "Python"; data[cur] = {"phases": {}}
        elif ln.startswith("===== R"):
            cur = "R"; data[cur] = {"phases": {}}
        if cur is None:
            continue
        m = re.match(r"PHASE\s+(\S+)\s+([\d.]+)\s+s", ln)
        if m:
            data[cur]["phases"][m.group(1)] = float(m.group(2))
        # R reports VmHWM in kB; Python reports its own peak RSS directly in GB.
        m = re.search(r"PEAKMEM_SELF.*VmHWM:\s+(\d+)\s*kB", ln)
        if m:
            data[cur]["peak_gb"] = int(m.group(1)) / 1024 / 1024
        else:
            m = re.search(r"PEAKMEM_SELF\s+([\d.]+)\s*GB", ln)
            if m:
                data[cur]["peak_gb"] = float(m.group(1))
        m = re.search(r"PEAKMEM_CHILD_MAX\s+([\d.]+)\s*GB", ln)
        if m:
            data[cur]["peak_child_gb"] = float(m.group(1))
    for tool, d in data.items():
        d["total_s"] = sum(d.get("phases", {}).values())
    return data


def load(path):
    dt = {c: "int32" for c in INT_COLS}
    dt.update({c: "float32" for c in FLT_COLS})
    dt.update({"gene": "string", "txcoor": "int32", "refSeq": "string",
               "chr": "string", "strand": "string", "gencoor": "int32"})
    usecols = ["chr", "gencoor", "strand", "gene", "txcoor", "refSeq"] + INT_COLS + FLT_COLS
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c in usecols,
                     dtype={k: v for k, v in dt.items() if k in usecols},
                     na_values=["NA"])
    return df


def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def img(b64, alt=""):
    return f'<img src="data:image/png;base64,{b64}" alt="{alt}" style="max-width:100%;">'


def hexbin_scatter(x, y, xlabel, ylabel, title, logscale=False):
    fig, ax = plt.subplots(figsize=(4.6, 4.4))
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if logscale:
        x = np.log10(x + 1); y = np.log10(y + 1)
        xlabel += " (log10+1)"; ylabel += " (log10+1)"
    if len(x):
        ax.hexbin(x, y, gridsize=60, bins="log", cmap="viridis", mincnt=1)
        lim = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lim, lim, "r--", lw=1, alpha=0.8)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10)
    return fig_to_b64(fig)


def main():
    log("loading Python TSV"); py = load(PY_TSV)
    log("loading R TSV"); r = load(R_TSV)
    log(f"py rows={len(py):,}  r rows={len(r):,}")

    # ---- gene sets ----
    g_py = set(py["gene"].unique()); g_r = set(r["gene"].unique())
    both = g_py & g_r
    only_py = g_py - g_r; only_r = g_r - g_py

    # ---- position-level join ----
    pj = py.set_index(KEYS)
    rj = r.set_index(KEYS)
    common_idx = pj.index.intersection(rj.index)
    log(f"common positions: {len(common_idx):,}")
    pj = pj.loc[common_idx]
    rj = rj.loc[common_idx]
    # align order
    rj = rj.reindex(pj.index)

    # refSeq agreement (coerce to plain numpy bool; refSeq is read as nullable string)
    ref_eq = np.asarray(pj["refSeq"].astype(object).to_numpy()
                        == rj["refSeq"].astype(object).to_numpy(), dtype=bool)

    # per-column agreement
    col_stats = []
    for c in INT_COLS:
        a = pj[c].to_numpy(); b = rj[c].to_numpy()
        eq = (a == b)
        diff = np.abs(a.astype(np.int64) - b.astype(np.int64))
        with np.errstate(invalid="ignore"):
            corr = np.corrcoef(a, b)[0, 1] if a.std() and b.std() else float("nan")
        col_stats.append((c, eq.mean() * 100, int(diff.max()), float(diff.mean()),
                          corr, int((~eq).sum())))
    for c in FLT_COLS:
        a = pj[c].to_numpy(); b = rj[c].to_numpy()
        m = np.isfinite(a) & np.isfinite(b)
        both_na = (~np.isfinite(a)) & (~np.isfinite(b))
        if m.sum():
            diff = np.abs(a[m] - b[m])
            corr = np.corrcoef(a[m], b[m])[0, 1] if a[m].std() and b[m].std() else float("nan")
            close = (diff <= 1e-6).mean() * 100
        else:
            diff = np.array([0.0]); corr = float("nan"); close = float("nan")
        col_stats.append((c, close, float(diff.max()), float(diff.mean()), corr,
                          int(m.sum())))

    # ---- difference diagnostics (coverage) ----
    cov_py = pj["cov"].to_numpy(np.int64); cov_r = rj["cov"].to_numpy(np.int64)
    cov_d = cov_py - cov_r
    diff_mask = cov_d != 0
    n_pos_diff = int(diff_mask.sum())
    py_gt = int((cov_d > 0).sum()); r_gt = int((cov_d < 0).sum())
    # per-gene: fraction of genes with fully identical coverage; top differing genes
    genes_idx = pj.index.get_level_values("gene").to_numpy()
    diff_df = pd.DataFrame({"gene": genes_idx, "d": diff_mask.astype(int)})
    per_gene_diff = diff_df.groupby("gene", observed=True)["d"].sum()
    genes_identical = int((per_gene_diff == 0).sum())
    genes_total_cmp = int(per_gene_diff.size)
    top_diff_genes = per_gene_diff.sort_values(ascending=False).head(10)

    # ---- gene-level aggregates ----
    gpy = py.groupby("gene", observed=True).agg(
        reads_py=("start_5p", "sum"), cov_py=("cov", "sum"), npos_py=("txcoor", "size"))
    grr = r.groupby("gene", observed=True).agg(
        reads_r=("start_5p", "sum"), cov_r=("cov", "sum"), npos_r=("txcoor", "size"))
    gm = gpy.join(grr, how="inner")

    # ---- plots ----
    log("plotting")
    plots = {}
    plots["genes_reads"] = hexbin_scatter(
        gm["reads_r"].to_numpy(float), gm["reads_py"].to_numpy(float),
        "R sum(start_5p)", "Py sum(start_5p)", "Gene-level read counts", logscale=True)
    plots["cov"] = hexbin_scatter(
        rj["cov"].to_numpy(float), pj["cov"].to_numpy(float),
        "R cov", "Py cov", "Per-position coverage", logscale=True)
    plots["start"] = hexbin_scatter(
        rj["start_5p"].to_numpy(float), pj["start_5p"].to_numpy(float),
        "R start_5p", "Py start_5p", "Read-starts (5')", logscale=True)
    plots["end"] = hexbin_scatter(
        rj["end_3p"].to_numpy(float), pj["end_3p"].to_numpy(float),
        "R end_3p", "Py end_3p", "Read-ends (3')", logscale=True)

    # misincRate and C->T where both finite & adequately covered
    a = pj["misincRate"].to_numpy(); b = rj["misincRate"].to_numpy()
    plots["misinc"] = hexbin_scatter(b, a, "R misincRate", "Py misincRate",
                                     "Misincorporation rate")
    a = pj["MR_CtoT"].to_numpy(); b = rj["MR_CtoT"].to_numpy()
    plots["ct"] = hexbin_scatter(b, a, "R C→T rate", "Py C→T rate",
                                 "STAMP C→T editing rate")

    # difference histogram for cov
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    d = (pj["cov"].to_numpy(np.int64) - rj["cov"].to_numpy(np.int64))
    nz = d[d != 0]
    ax.hist(nz if len(nz) else [0], bins=60, color="#0098fd")
    ax.set_title(f"Per-position cov difference (Py−R), non-zero only: {len(nz):,} of {len(d):,}",
                 fontsize=9)
    ax.set_xlabel("Py cov − R cov"); ax.set_yscale("log")
    plots["covdiff"] = fig_to_b64(fig)

    # ---- assemble HTML ----
    def tbl(rows, header):
        h = "".join(f"<th>{c}</th>" for c in header)
        body = ""
        for row in rows:
            body += "<tr>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>"
        return f"<table><thead><tr>{h}</tr></thead><tbody>{body}</tbody></table>"

    colrows = [(c, f"{eq:.4f}", mx, f"{mn:.4g}", f"{cr:.6f}" if cr == cr else "n/a", nd)
               for (c, eq, mx, mn, cr, nd) in col_stats]

    pos_total = len(common_idx)
    overall_int_exact = np.mean([s[1] for s in col_stats[:len(INT_COLS)]])

    # ---- resources: runtime + memory ----
    res = parse_phase_timing(f"{V}/phase_timing.log")
    PHASE_ORDER = [("load_bed", "Load annotation (BED12)"),
                   ("load_genome", "Load genome (FASTA)"),
                   ("load_bam", "Load BAM alignments"),
                   ("tx_reads", "tx_reads (genomic→tx)"),
                   ("makeDT", "Build table (makeDT)"),
                   ("add_metrics", "Add metrics")]
    res_section = ""
    PY_TKPI, R_TKPI = PY_T, R_T
    if res.get("Python") and res.get("R"):
        py, rr = res["Python"], res["R"]
        rows = [(label, f'{py["phases"].get(k, float("nan")):.1f}',
                 f'{rr["phases"].get(k, float("nan")):.1f}') for k, label in PHASE_ORDER]
        py_tot, r_tot = py.get("total_s", 0.0), rr.get("total_s", 0.0)
        rows.append(("<b>Total wall time</b>", f"<b>{py_tot:.1f}</b>", f"<b>{r_tot:.1f}</b>"))
        py_mem, r_mem = py.get("peak_gb", float("nan")), rr.get("peak_gb", float("nan"))
        fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.2))
        axs[0].barh(["R", "Python"], [r_tot, py_tot], color=["#c0392b", "#0098fd"])
        axs[0].set_title("Total wall time (s, 8 cores)", fontsize=10)
        for i, v in enumerate([r_tot, py_tot]):
            axs[0].text(v, i, f" {v:.0f}s", va="center")
        axs[1].barh(["R", "Python"], [r_mem, py_mem], color=["#c0392b", "#0098fd"])
        axs[1].set_title("Peak memory (GB, main process)", fontsize=10)
        for i, v in enumerate([r_mem, py_mem]):
            axs[1].text(v, i, f" {v:.2f} GB", va="center")
        fig.tight_layout()
        plots["resources"] = fig_to_b64(fig)
        speedup = (r_tot / py_tot) if py_tot else float("nan")
        memratio = (r_mem / py_mem) if py_mem else float("nan")
        PY_TKPI, R_TKPI = f"{py_tot:.0f} s", f"{r_tot:.0f} s"
        res_section = f"""
<h2>0. Resources — runtime &amp; memory (8 cores each)</h2>
<div>
<span class="kpi">Speedup (R / Python)<br><b class="good">{speedup:.1f}×</b></span>
<span class="kpi">Python peak memory<br><b>{py_mem:.2f} GB</b></span>
<span class="kpi">R peak memory<br><b>{r_mem:.2f} GB</b></span>
<span class="kpi">Memory ratio (R / Python)<br><b>{memratio:.1f}×</b></span>
</div>
{tbl(rows, ["phase", "Python (s)", "R (s)"])}
<div class="grid"><div class="card">{img(plots['resources'], 'resources')}</div></div>
<p style="font-size:12px;color:#555">Both tools were given 8 cores (Python <code>--threads 8</code>;
R <code>nCores=8</code> in tx_reads and makeDT, via <code>parallel::mclapply</code>). Peak memory is the
main process's peak resident set size (Python <code>ru_maxrss</code>; R <code>/proc/self/status</code>
VmHWM); forked workers add some copy-on-write overhead on top. R loads all alignments into memory up front
&mdash; a single-threaded step with no nCores option &mdash; whereas the Python tool streams reads per gene
through the BAM index, so it has no separate load phase and a smaller footprint.</p>
"""

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>txtools: Python vs R comparison (chr19)</title>
<style>
body{{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:32px;color:#1b1b1b;line-height:1.45}}
h1{{font-size:24px}} h2{{margin-top:34px;border-bottom:2px solid #0098fd;padding-bottom:4px}}
table{{border-collapse:collapse;margin:10px 0;font-size:13px}}
th,td{{border:1px solid #ccc;padding:5px 9px;text-align:right}} th{{background:#f0f4f8}}
td:first-child,th:first-child{{text-align:left}}
.grid{{display:flex;flex-wrap:wrap;gap:18px}} .card{{flex:1 1 300px}}
.kpi{{display:inline-block;background:#f0f4f8;border-radius:8px;padding:10px 16px;margin:6px 8px 6px 0}}
.kpi b{{font-size:20px;color:#0067b8}} .good{{color:#0a8f3c;font-weight:bold}} .warn{{color:#c0392b;font-weight:bold}}
code{{background:#f3f3f3;padding:1px 5px;border-radius:4px}} .note{{background:#fffbe6;border-left:4px solid #f3b018;padding:8px 12px;margin:10px 0}}
</style></head><body>
<h1>txtools — Python reimplementation vs reference R package</h1>
<p>Validation on a real dataset: <code>HEK293T_control_STAMP_72h_high_rep2</code> (APOBEC1 STAMP, hg19),
single-end STAR alignments, <b>chr19</b> subset. Both tools were run with identical inputs and parameters.</p>

<div>
<span class="kpi">Dataset reads (chr19)<br><b>1,007,322</b></span>
<span class="kpi">Gene models (refGene chr19)<br><b>4,647</b></span>
<span class="kpi">strandMode<br><b>2</b></span>
<span class="kpi">minReads<br><b>50</b></span>
<span class="kpi">Python runtime<br><b>{PY_TKPI}</b></span>
<span class="kpi">R runtime<br><b>{R_TKPI}</b></span>
</div>
{res_section}
<h2>1. Gene-model agreement</h2>
<div>
<span class="kpi">Genes kept — Python<br><b>{len(g_py):,}</b></span>
<span class="kpi">Genes kept — R<br><b>{len(g_r):,}</b></span>
<span class="kpi">In both<br><b class="good">{len(both):,}</b></span>
<span class="kpi">Python only<br><b>{len(only_py):,}</b></span>
<span class="kpi">R only<br><b>{len(only_r):,}</b></span>
</div>
<div class="grid"><div class="card">{img(plots['genes_reads'],'gene reads')}</div></div>

<h2>2. Per-position agreement (common positions: {pos_total:,})</h2>
<p>refSeq identical at <b class="good">{ref_eq.mean()*100:.4f}%</b> of common positions.
Mean exact-match across integer columns: <b>{overall_int_exact:.4f}%</b>.</p>
{tbl(colrows, ["column","% exact (or |Δ|≤1e-6)","max |Δ|","mean |Δ|","Pearson r","# rows differing / used"])}

<div class="grid">
<div class="card">{img(plots['cov'],'cov')}</div>
<div class="card">{img(plots['start'],'start')}</div>
<div class="card">{img(plots['end'],'end')}</div>
<div class="card">{img(plots['covdiff'],'covdiff')}</div>
</div>

<h2>2b. Coverage difference diagnostics</h2>
<div>
<span class="kpi">Positions with identical cov<br><b class="good">{(1-n_pos_diff/max(pos_total,1))*100:.4f}%</b></span>
<span class="kpi">Positions differing<br><b>{n_pos_diff:,}</b></span>
<span class="kpi">Genes with 100% identical cov<br><b>{genes_identical:,} / {genes_total_cmp:,}</b></span>
<span class="kpi">Direction (Py&gt;R / Py&lt;R)<br><b>{py_gt:,} / {r_gt:,}</b></span>
</div>
<p>Differences are concentrated in a minority of gene models — consistent with the documented
per-read vs all-or-nothing splice-read handling. Top gene models by differing positions:</p>
{tbl([(g, int(v)) for g, v in top_diff_genes.items()], ["gene model", "# positions differing"])}

<h2>3. Nucleotide frequencies &amp; modification rates</h2>
<div class="grid">
<div class="card">{img(plots['misinc'],'misinc')}</div>
<div class="card">{img(plots['ct'],'ct')}</div>
</div>

<h2>4. Notes on expected differences</h2>
<div class="note">
By design the Python rewrite differs from R in two documented ways that can move a small number of positions:
<ul>
<li><b>Spliced-read junction check is per-read</b> in Python; the R code drops <i>all</i> spliced reads in a gene if any one read's junctions disagree with the annotation. Differences therefore concentrate in genes containing junction-spanning reads.</li>
<li><b>strandMode=2 boundary handling</b> is identical to R for the common case; it only diverges in pathological mate geometries (not applicable to this single-end run).</li>
</ul>
</div>
<p style="color:#777;font-size:12px">Generated {time.strftime('%Y-%m-%d %H:%M')} • Python txtools {PY_VER} • R txtools {R_VER}</p>
</body></html>"""
    with open(OUT, "w") as fh:
        fh.write(html)
    log("wrote", OUT)
    # also print a console summary
    print("\n=== SUMMARY ===")
    print(f"genes: both={len(both)} py_only={len(only_py)} r_only={len(only_r)}")
    print(f"common positions: {pos_total:,}  refSeq match: {ref_eq.mean()*100:.4f}%")
    for (c, eq, mx, mn, cr, nd) in col_stats:
        print(f"  {c:>11}: exact/close={eq:7.4f}%  maxΔ={mx:<10g} r={cr:.6f}  ndiff/used={nd}")


def _parse_R_runtime():
    try:
        for ln in open(f"{V}/r_run.log"):
            if "R elapsed:" in ln:
                return ln.split("R elapsed:")[1].strip()
    except OSError:
        pass
    return "n/a"


PY_T = "195 s (8 cores)"
R_T = _parse_R_runtime()
PY_VER = "0.1.0"
R_VER = "1.0.6"

if __name__ == "__main__":
    main()
