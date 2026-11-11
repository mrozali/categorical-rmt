"""S3 qutrit chain, Sec. 2.4.3: Table 1 and Figs. 3, 5 and 6 from the two stored ensembles (grouped wrapper).

Runs code/s3/python/run_paper_numerics.py, unchanged, on
  data/s3/S3_rep_category_L8_spectra.npz  (June 2026 ensemble: Table 1, Fig. 3)  -> generated/s3_june/
  data/s3/ensemble_L8_n120.npz            (November 2026 ensemble: Figs. 5, 6)   -> generated/s3_nov/
then checks the outputs against the paper:
  * Table 1: every <r> and its standard error, rounded to 4 decimals, equal the values parsed from
    paper/categorical_rmt.tex (June ensemble);
  * block dimensions 1093 to 3281, and 3^8 = 6561 per twisted Hilbert space;
  * the Theta pair ([c],omega), ([c],omega^2): equal spectra (< 1e-11) and C_AB = 1 (> 0.9999);
  * Fig. 5 and Fig. 6 (November ensemble): pixel comparison with the paper's PNGs, PASS up to anti-aliasing
    (max channel difference <= 0.02);
  * Fig. 3 is written for visual comparison only. Its legend gap counts differ from the paper by up to 8,
    because the paper's version (June plotting script, not included) cut gaps at s < 5 after its own unfolding.
It also prints, without PASS criteria, the offsets of Table 1 from the RMT benchmarks (footnote of Sec. 2.4.3),
the ratio K_c/K_RMT in bands of t, and the off-diagonal C_AB. Fig. 4 has its own script,
fig04_s3_sff_goe_binned.py. About 10 s. Last line: RESULT: PASS or RESULT: FAIL ...
"""
import json, os, re, subprocess, sys
import numpy as np
from matplotlib.image import imread

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
GEN = os.path.join(ROOT, "code", "figure-reproduction", "generated")
RUN = os.path.join(ROOT, "code", "s3", "python", "run_paper_numerics.py")
RUNS = {"june": os.path.join(ROOT, "data", "s3", "S3_rep_category_L8_spectra.npz"),
        "nov": os.path.join(ROOT, "data", "s3", "ensemble_L8_n120.npz")}
SECTORS = ["([e],1)", "([e],sgn)", "([e],2)", "([t],+)", "([t],-)", "([c],1)", "([c],omega)", "([c],omega^2)"]
TEX_NAME = {"([e],1)": r"([e],1)", "([e],sgn)": r"([e],\sgn)", "([e],2)": r"([e],2)", "([t],+)": r"([t],+)",
            "([t],-)": r"([t],-)", "([c],1)": r"([c],1)", "([c],omega)": r"([c],\omega)", "([c],omega^2)": r"([c],\omega^2)"}
BENCH = {"GOE": 0.5307, "GUE": 0.5996}
fails = []

# Table 1 as printed in the paper
tex = open(os.path.join(ROOT, "paper", "categorical_rmt.tex"), encoding="utf-8").read()
tab = tex[tex.index(r"\label{tab:S3r}") - 1500: tex.index(r"\label{tab:S3r}")]
paper = {}
for s in SECTORS:
    m = re.search(r"\$" + re.escape(TEX_NAME[s]) + r"\$\s*&[^&]*&\s*\$([\d.]+)\s*\\pm\s*([\d.]+)\$\s*&\s*(\w+)", tab)
    if not m:
        sys.exit(f"RESULT: FAIL could not parse Table 1 row {s} from the paper")
    paper[s] = (float(m[1]), float(m[2]), m[3])

out = {}
for name, npz in RUNS.items():
    d = os.path.join(GEN, f"s3_{name}")
    p = subprocess.run([sys.executable, "-B", RUN, "--from-npz", npz, "--out", d], capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"RESULT: FAIL run_paper_numerics.py on {os.path.basename(npz)}:\n{p.stdout}\n{p.stderr}")
    out[name] = json.load(open(os.path.join(d, "numerics_summary.json")))
    print(f"ran run_paper_numerics.py on data/s3/{os.path.basename(npz)} -> generated/s3_{name}/")

# Table 1
t1 = out["june"]["table1"]
print("\nTable 1 (June ensemble), <r> +- SEM: this run | paper | November ensemble, for comparison")
for s in SECTORS:
    r, e = t1[s]["rbar"], t1[s]["sem"]
    pr, pe, cls = paper[s]
    ok = round(r, 4) == pr and round(e, 4) == pe and t1[s]["class"] == cls
    if not ok:
        fails.append(f"Table 1 {s}: {r:.4f}+-{e:.4f} vs paper {pr}+-{pe} {cls}")
    n = out["nov"]["table1"][s]
    print(f"  {s:14s} {r:.4f} +- {e:.4f} | {pr:.4f} +- {pe:.4f} {cls} {'ok' if ok else 'DIFFERS'} | {n['rbar']:.4f} +- {n['sem']:.4f}")

z = [abs(t1[s]["rbar"] - out["nov"]["table1"][s]["rbar"]) / np.hypot(t1[s]["sem"], out["nov"]["table1"][s]["sem"])
     for s in SECTORS]
print(f"  November against June: largest difference {max(z):.1f} combined standard errors")

dims = out["june"]["dimensions"]; per = out["june"]["twisted_hilbert_space_dims"]
print(f"\nblock dimensions {min(dims.values())} to {max(dims.values())}, twisted Hilbert spaces {sorted(set(per.values()))}"
      "   (paper: 1093 to 3281, 3^8 = 6561)")
if (min(dims.values()), max(dims.values())) != (1093, 3281) or set(per.values()) != {6561}:
    fails.append("dimensions")

full = [abs(t1[s]["rbar"] - BENCH[t1[s]["class"]]) / t1[s]["sem"] for s in SECTORS]
rnd = [abs(paper[s][0] - BENCH[paper[s][2]]) / paper[s][1] for s in SECTORS]
print(f"offset below the RMT benchmark, in standard errors: {min(full):.1f} to {max(full):.1f} from the unrounded values, "
      f"{min(rnd):.1f} to {max(rnd):.1f} from the rounded Table 1 entries   (paper footnote: 1.5 to 5.3, arXiv v1: 1.5 to 5.7)")

# the Theta pair and C_AB
for name in RUNS:
    o = out[name]; C = o["C"]
    print(f"\n{name}: Theta pair max |dE| = {o['theta_max_deviation']:.2e}, cross vs auto SFF rms on 0.2<=t<=3 = "
          f"{o['cross_vs_auto_rms_0.2to3']:.1e}, C_AB(Theta pair) = {C['theta_pair']:.4f}, "
          f"other 27 pairs max |C_AB| = {C['max_off_diagonal']:.3f} (rms {C['rms_off_diagonal']:.3f})")
    if o["theta_max_deviation"] >= 1e-11 or C["theta_pair"] <= 0.9999:
        fails.append(f"{name}: Theta pair")
    bands = o["sff_ratio_bands"]
    for b in ("0.02-0.05", "0.2-1.0"):
        v = [bands[s][b] for s in SECTORS]
        print(f"  K_c/K_RMT averaged over {b} t_H: {min(v):.2f} to {max(v):.2f}"
              + ("   (paper: RMT for t >~ 0.2 t_H)" if b == "0.2-1.0" else "   (paper: above the ramp at early times)"))

# figures
def compare(gen, ref):
    A, B = imread(gen)[..., :3], imread(ref)[..., :3]
    if A.shape != B.shape:
        return None, f"shape {A.shape} vs paper {B.shape}"
    dd = np.abs(A - B).max(axis=2)
    return float(dd.max()), f"max channel difference {dd.max():.3f}, differing pixels {(dd > 0).mean():.5f}"
print()
for fig, gen, ref in (("Fig. 5", "s3_nov/fig5_sff_GUE.png", "S3_rep_category_L8_sff_GUE.png"),
                      ("Fig. 6", "s3_nov/fig6_correlation_matrix.png", "S3_rep_category_L8_cross_sector.png")):
    mx, msg = compare(os.path.join(GEN, gen), os.path.join(ROOT, "paper", ref))
    print(f"{fig}: generated/{gen} vs paper/{ref}: {msg}")
    if mx is None or mx > 0.02:
        fails.append(f"{fig} pixels")
print("Fig. 3: generated/s3_june/fig3_spacing.png, compare visually with paper/S3_rep_category_L8_spacing_hist.png")
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
