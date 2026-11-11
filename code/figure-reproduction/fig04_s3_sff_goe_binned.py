"""Fig. 4 (fig:sffGOE): connected SFF of the six GOE-class S3 sectors, on the paper's time grid.

The paper's Fig. 4 was drawn by the authors' June 2026 plotting script, which is not included. According to the
documentation of code/s3/python/, it used 400 log-spaced times from 10^-1.7 to 10^0.5 averaged in 44 logarithmic
bins. This script redraws it from the stored June ensemble with the tools of code/s3/python/, on that grid.
The binning convention is reconstructed: K_c is averaged over the times in each bin and plotted at their geometric
mean. The result has the same shape as the paper's figure but does not match it point for point, so the status is
runs-but-differs. code/s3/python/run_paper_numerics.py draws an unbinned 80-time version.

Writes generated/fig04_s3_sff_GOE_binned.png and prints K_c/K_GOE averaged over 0.2 <= t <= 1 for each sector
(caption: "matching the GOE ramp -> plateau for t >~ 0.2 t_H") and over 0.02 <= t <= 0.05.
PASS if every 0.2-1 average lies within 10% of 1. About 10 s. Last line: RESULT: PASS or RESULT: FAIL ...
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
GEN = os.path.join(ROOT, "code", "figure-reproduction", "generated")
sys.path.insert(0, os.path.join(ROOT, "code", "s3", "python"))
sys.dont_write_bytecode = True
from run_paper_numerics import load_npz, LABEL, PLAIN
from s3_chain import SECTORS, SECTOR_CLASS
from spectral_tools import amplitudes, sff_connected, sff_rmt

spec, _, _, L = load_npz(os.path.join(ROOT, "data", "s3", "S3_rep_category_L8_spectra.npz"))
ts = np.logspace(-1.7, 0.5, 400)
edges = np.logspace(-1.7, 0.5, 45)
idx = np.clip(np.digitize(ts, edges) - 1, 0, 43)
fig, ax = plt.subplots(figsize=(8.0, 5.6))
fails = []
print("K_c/K_GOE averaged over the 400-time grid:   0.2<=t<=1   0.02<=t<=0.05")
for k, s in enumerate([s for s in SECTORS if SECTOR_CLASS[s] == "GOE"]):
    Z, Ns = amplitudes(spec[s], ts)
    K = sff_connected(Z, Ns)
    ratio = K / sff_rmt(ts, "GOE")
    late, early = ratio[(ts >= 0.2) & (ts <= 1)].mean(), ratio[(ts >= 0.02) & (ts <= 0.05)].mean()
    print(f"  {PLAIN[s]:14s} dim {spec[s].shape[1]:5d}   {late:6.3f}      {early:6.2f}")
    if abs(late - 1) > 0.1:
        fails.append(f"{PLAIN[s]} 0.2-1 average {late:.3f}")
    bins = [b for b in range(44) if np.any(idx == b)]
    tb = np.array([np.exp(np.log(ts[idx == b]).mean()) for b in bins])
    Kb = np.array([K[idx == b].mean() for b in bins])
    ax.loglog(tb, Kb, "o-", color=f"C{k}", ms=3.5, lw=1.0, label=LABEL[s])
tg = np.logspace(np.log10(0.02), np.log10(3.0), 200)
ax.loglog(tg, sff_rmt(tg, "GOE"), "k-", lw=1.8, label="GOE theory")
ax.loglog(tg, sff_rmt(tg, "GUE"), "r--", lw=1.2, label="GUE theory")
ax.axvline(1.0, color="gray", ls=":", lw=1); ax.set_xlabel("t  (Heisenberg time = 1)")
ax.set_ylabel(r"connected SFF  $K_c(t)/N$")
ax.set_title(f"Spectral form factor of the 6 GOE-class sectors (L = {L}, {spec['e1'].shape[0]} realizations)", fontsize=11)
ax.legend(fontsize=9, loc="upper left", ncol=2); fig.tight_layout()
os.makedirs(GEN, exist_ok=True)
out = os.path.join(GEN, "fig04_s3_sff_GOE_binned.png")
fig.savefig(out, dpi=150); plt.close(fig)
print(f"wrote generated/{os.path.basename(out)}; compare visually with paper/S3_rep_category_L8_sff_GOE.png")
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
