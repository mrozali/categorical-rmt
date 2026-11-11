"""Fig. 7 reproduction: L=26 Fibonacci chain, lambda=1.25, polynomial-unfolded P(s) in the
five naive blocks (untwisted W+, W-, twisted (tau,tau), (tau,1), (1,tau)).

Reads  <repo root>/data/fib/fib_L26_lam1.25_spectra.npz
Writes <repo root>/code/figure-reproduction/generated/fib_sectors_L26.{pdf,png}
Prints the five <r> values with SEM, checks they round to 0.5292, 0.5309, 0.5309, 0.5306, 0.5306,
then (if pdftoppm is available) compares a 50 dpi render of the generated PDF with
paper/fib_sectors_L26.pdf. Last line: RESULT: PASS or RESULT: FAIL ...

Run from anywhere:  python code/figure-reproduction/fig07_fib_sectors_L26.py
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SPEC_FILE = os.path.join(ROOT, "data", "fib", "fib_L26_lam1.25_spectra.npz")
FIG_DIR = os.path.join(HERE, "generated")

UNTW_SECTORS = ["W+", "W-"]
TW_SECTORS = ["W1", "W2", "W3"]
ALL_SECTORS = UNTW_SECTORS + TW_SECTORS

L_MAIN = 26
MIN_LEVELS = 13

COLORS = {
    "W+": "seagreen", "W-": "darkorange",
    "W1": "mediumpurple", "W2": "cadetblue", "W3": "rosybrown",
}

# Panel titles: untwisted pieces W_+/W_- (eigenvalues +phi, -1/phi of L^{tau,tau}_{1<-1});
# twisted pieces by center label. Twisted sector s<->m mod 5 (eth_twisted.py): 0 -> (tau,tau),
# 2 -> spin e^{4 pi i/5} = (tau,1), 3 -> spin e^{-4 pi i/5} = (1,tau)  [paper eq:spins].
TITLES = {
    "W+": r"untwisted, $W_+$",
    "W-": r"untwisted, $W_-$",
    "W1": r"$\tau$-twisted, $(\tau,\tau)$",
    "W2": r"$\tau$-twisted, $(\tau,1)$",
    "W3": r"$\tau$-twisted, $(1,\tau)$",
}


# -- verbatim from the original (private) fib_unfold_replot.py
# (duplicated from catrmt_fig_fib_sectors.py, same directory) --
def unfold_local(e, deg=7, edge_drop=0.05):
    e = np.sort(np.real(e)); n = len(e)
    if n < deg + 5: return np.array([])
    lo = int(edge_drop * n); hi = n - lo
    e_in = e[lo:hi]; k_in = np.arange(lo, hi, dtype=float)
    coeffs = np.polyfit(e_in, k_in, deg)
    xi = np.polyval(coeffs, e_in)
    s = np.diff(xi); s = s[s > 0]
    if len(s) == 0: return s
    return s / np.mean(s)

def ratios_local(e, edge_drop=0.05):
    e = np.sort(np.real(e)); n = len(e)
    lo = int(edge_drop * n); hi = n - lo
    if hi - lo < 6: return np.array([])
    s = np.diff(e[lo:hi])
    if len(s) < 2: return np.array([])
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r[np.isfinite(r)]
# -- end verbatim block --


# -- data loaders and pooling, duplicated from catrmt_fig_fib_sectors.py --
def _lean_untw_E(L):
    """Adapted: read the compact file (keys 'U|<orig key>') instead of the Fir bank; same
    self-conjugate momentum hygiene as eth_fir_load.fir_flat_untw."""
    assert L == 26
    out = {}
    with np.load(SPEC_FILE) as z:
        for kk in z.files:
            if not kk.startswith("U|"):
                continue
            k = kk[2:]
            m = re.match(r"k(\d+)W", k)
            if m and int(m.group(1)) in (0, L // 2):
                continue
            out[k] = z[kk]
    return out or None


def _lean_tw_E(L):
    """Adapted: read the compact file (keys 'T|<orig key>'), same sector hygiene (a sector with
    fewer than L/2 blocks is dropped)."""
    assert L == 26
    out = {}
    with np.load(SPEC_FILE) as z:
        for kk in z.files:
            if kk.startswith("T|"):
                out[kk[2:]] = z[kk]
    labs = {k.split("|")[0] for k in out}
    counts = {s: sum(1 for l in labs if l.startswith(f"s{s}")) for s in "123"}
    keep = {s for s, c in counts.items() if c >= L // 2}
    if not keep:
        return None
    return {k: v for k, v in out.items() if k.split("|")[0][1] in keep}


def untw_blocks(L, sector):
    """List of eigenvalue arrays for untwisted sector ('W+'/'W-') at size L."""
    sign = re.escape(sector[-1])
    regex = re.compile(r"k\d+W%s\|E" % sign)
    d = _lean_untw_E(L)
    if not d:
        return []
    return [v for k, v in d.items() if regex.fullmatch(k)]


def tw_blocks(L, sector):
    """List of eigenvalue arrays for twisted sector ('W1'/'W2'/'W3') at size L,
    after sector-1 monodromy-block hygiene (drop m=0, m>=5L/2)."""
    sidx = sector[-1]
    regex = re.compile(r"s(\d)m(\d+)\|E")
    d = _lean_tw_E(L)
    if not d:
        return []
    items = [(k, v) for k, v in d.items() if regex.fullmatch(k)]
    mmax = 5 * L // 2
    out = []
    for k, e in items:
        mm = regex.fullmatch(k)
        s, m = mm.group(1), int(mm.group(2))
        if s != sidx:
            continue
        if sidx == "1" and (m == 0 or m >= mmax):
            continue
        out.append(e)
    return out


def sector_blocks(L, sector):
    if sector in UNTW_SECTORS:
        return untw_blocks(L, sector)
    return tw_blocks(L, sector)


def pooled_stats(blocks_e):
    """Pool unfolded spacings and gap ratios across blocks (>=13 levels each).
    Returns (spacings, ratios, n_blocks_used)."""
    spacings, ratios, n_used = [], [], 0
    for e in blocks_e:
        if len(e) < MIN_LEVELS:
            continue
        n_used += 1
        s = unfold_local(e, deg=7, edge_drop=0.05)
        if len(s):
            spacings.append(s)
        r = ratios_local(e)
        if len(r):
            ratios.append(r)
    spacings = np.concatenate(spacings) if spacings else np.array([])
    ratios = np.concatenate(ratios) if ratios else np.array([])
    return spacings, ratios, n_used


def mean_sem(x):
    x = np.asarray(x)
    if len(x) == 0:
        return np.nan, np.nan
    return float(np.mean(x)), float(np.std(x) / np.sqrt(len(x)))
# -- end duplicated block --


# ---------------------------------------------------------------------------
# Gather L=26 data for the five P(s) panels
# ---------------------------------------------------------------------------
main_data = {}
print(f"=== L={L_MAIN} sector panels ===")
for sector in ALL_SECTORS:
    blocks = sector_blocks(L_MAIN, sector)
    spacings, ratios, n_used = pooled_stats(blocks)
    rmean, rsem = mean_sem(ratios)
    main_data[sector] = dict(spacings=spacings, n_spacings=len(spacings),
                              rmean=rmean, rsem=rsem, n_used=n_used)
    print(f"{sector}: n_blocks_used={n_used}, n_spacings={len(spacings)}, "
          f"mean_r={rmean:.4f}, sem_r={rsem:.4f}")

for sector in ALL_SECTORS:
    rmean = main_data[sector]["rmean"]
    assert 0.52 <= rmean <= 0.54, (
        f"L={L_MAIN} sector {sector} mean gap ratio {rmean:.4f} outside "
        f"expected GOE window [0.52, 0.54] -- aborting."
    )

# ---------------------------------------------------------------------------
# Figure (matplotlib defaults, paper style copied from make_L22_clean_plot.py)
# ---------------------------------------------------------------------------
s_grid = np.linspace(0.001, 4, 500)
GOE = (np.pi / 2) * s_grid * np.exp(-np.pi / 4 * s_grid ** 2)
Pois = np.exp(-s_grid)

fig = plt.figure(figsize=(14, 9.5))
gs = fig.add_gridspec(2, 6)

PANEL_SLOTS = {
    "W+": gs[0, 1:3], "W-": gs[0, 3:5],
    "W1": gs[1, 0:2], "W2": gs[1, 2:4], "W3": gs[1, 4:6],
}

for i, sector in enumerate(ALL_SECTORS):
    ax = fig.add_subplot(PANEL_SLOTS[sector])
    d = main_data[sector]
    ax.hist(d["spacings"], bins=80, range=(0, 4), density=True, alpha=0.6,
            color=COLORS[sector],
            label=f"{d['n_spacings']} spacings,  $\\langle r\\rangle={d['rmean']:.4f}$")
    goe_label = r"single GOE  ($\langle r\rangle=0.5307$)" if i == 0 else "single GOE"
    ax.plot(s_grid, GOE, "b-", lw=2.5, label=goe_label)
    ax.plot(s_grid, Pois, "r--", lw=1.4, label="Poisson")
    ax.set_xlim(0, 4); ax.set_ylim(0, 0.95)
    ax.set_xlabel("s", fontsize=12); ax.set_ylabel("P(s)", fontsize=12)
    ax.set_title(TITLES[sector], fontsize=14)
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(alpha=0.25)

plt.suptitle(
    r"L=26 Fibonacci chain, $\lambda$=1.25:  polynomial-unfolded $P(s)$ in the "
    r"five naive blocks",
    fontsize=13)
plt.tight_layout(rect=[0, 0, 1, 0.96])

pdf_path = os.path.join(FIG_DIR, "fib_sectors_L26.pdf")
png_path = os.path.join(FIG_DIR, "fib_sectors_L26.png")
os.makedirs(FIG_DIR, exist_ok=True)
plt.savefig(pdf_path, bbox_inches="tight")
plt.savefig(png_path, dpi=140, bbox_inches="tight")

print("\nSaved:")
print(os.path.abspath(pdf_path))
print(os.path.abspath(png_path))


# ---------------------------------------------------------------------------
# Checks (added for the reproduction package)
# ---------------------------------------------------------------------------
EXPECTED = {"W+": "0.5292", "W-": "0.5309", "W1": "0.5309", "W2": "0.5306", "W3": "0.5306"}
print("\n<r> (pooled gap ratio, mean +- SEM):")
fails = []
for sector in ALL_SECTORS:
    d = main_data[sector]
    got = f"{d['rmean']:.4f}"
    ok = got == EXPECTED[sector]
    print(f"  {sector}: {d['rmean']:.6f} +- {d['rsem']:.4f}  rounds to {got}, expected {EXPECTED[sector]}  {'ok' if ok else 'MISMATCH'}")
    if not ok:
        fails.append(f"<r> {sector} {got} != {EXPECTED[sector]}")

paper_pdf = os.path.join(ROOT, "paper", "fib_sectors_L26.pdf")
PDFTOPPM = shutil.which("pdftoppm") or "/opt/homebrew/bin/pdftoppm"
pix = "skipped (pdftoppm not found)"
if os.path.exists(PDFTOPPM) and os.path.exists(paper_pdf):
    from matplotlib.image import imread
    with tempfile.TemporaryDirectory() as td:
        for tag, pdf in (("gen", pdf_path), ("paper", paper_pdf)):
            subprocess.run([PDFTOPPM, "-r", "50", "-png", pdf, os.path.join(td, tag)], check=True)
        def _load(tag):
            fn = [f for f in os.listdir(td) if f.startswith(tag)][0]
            return imread(os.path.join(td, fn))[..., :3]
        a, b = _load("gen"), _load("paper")
        if a.shape != b.shape:
            pix = f"shape mismatch {a.shape} vs {b.shape}"
            fails.append("pixel check " + pix)
        else:
            diff = np.abs(a.astype(float) - b.astype(float))
            frac = float(np.mean(diff.max(axis=-1) > 0))
            pix = f"max diff {diff.max():.4f}, fraction of differing pixels {frac:.5f}, shape {a.shape}"
            if diff.max() > 0.05 or frac > 0.01:
                fails.append("pixel check " + pix)
print("pixel check (50 dpi, generated vs paper/fib_sectors_L26.pdf):", pix)
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
