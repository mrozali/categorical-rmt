"""
Reproduces every number and figure of the S3 test in Section 2.4.3 of "Categorical RMT": L = 8, representative
twists e, sigma, r, 120 disorder realizations, central 70% of every spectrum, degree-8 unfolding.  The form
factors are evaluated on 80 log-spaced times from 0.02 to 3 Heisenberg times and C_AB is averaged over
0.1 <= t <= 1; these are the values used to make the manuscript's figures (the text itself only says "a band
of unfolded times").

Usage (the working directory does not matter; see DOCUMENTATION.md for details)
    python3 run_paper_numerics.py                       full run, output in results/
    python3 run_paper_numerics.py --seed-base 20248601  the June ensemble behind Table 1 and Figures 3-4
    python3 run_paper_numerics.py --from-npz FILE       analyse a stored ensemble (format of this script, or the
                                                        June file S3_rep_category_L8_spectra.npz)
Options (the window and the unfolding degree are the manuscript's; the other defaults are those of the run
behind Figures 5 and 6)
    --nreal N            number of disorder realizations, at least 2 (120)
    --seed-base S        realization r = 0, ..., N-1 uses numpy default_rng(S + r) (20261103)
    --L L                chain length, at least 4 (8); ignored with --from-npz, which takes L from the file
    --fraction F         central fraction of every spectrum used by all diagnostics (0.7)
    --degree D           degree of the polynomial that unfolds the spectra (8)
    --times TMIN TMAX N  N log-spaced times, in Heisenberg times, for the form factors (0.02 3 80)
    --band TMIN TMAX     time band over which C_AB is averaged, containing at least one grid time (0.1 1)
    --out DIR            output folder
A value out of range stops the script with an error message before any work.
Output folder: results/ for the default full run, otherwise results_<non-default settings>/ (for example
results_n10/ or results_n120_seed20248601/).  --from-npz FILE with the default analysis writes next to FILE,
replacing the tables and figures there, and with other settings to results_<file name>_<tags>/.

The default seed base 20261103 is the one of the ensemble behind Figures 5 and 6 as they now appear in the
manuscript, so the default run reproduces those two figures (see README.md for the check).  Table 1 and
Figures 3 and 4 come from the June 2026 ensemble, whose seed rule 20240601 + 1000 L + r is seed base 20248601
at L = 8: --seed-base 20248601 regenerates it and reproduces Table 1 exactly.

Output files (in results/ for the full run):
    run_log.txt                                  everything printed by this script (redraw_log.txt for --from-npz)
    ensemble_L8_n<nreal>.npz                     all sector spectra and couplings (16 MB for 120 realizations)
    table1_gap_ratios.json, table1_gap_ratios.png   Table 1 (<r> per sector, raw gaps) and its plot
    fig3_spacing.png                             Figure 3 (P(s) per sector)
    fig4_sff_GOE.png                             Figure 4 (connected SFF, six GOE-class sectors)
    fig5_sff_GUE.png                             Figure 5 (SFF and cross SFF of the Theta pair)
    fig6_correlation_matrix.png                  Figure 6 (C_AB)
    numerics_summary.json                        the main results and the settings used
Measured on an Apple M4 (macOS 15.6, Python 3.9.6, numpy 2.0.2): 1057 s for the full run, 4 s for --from-npz.
Definitions: s3_chain.py (model, eqs. (2.23)-(2.24), sector bases) and spectral_tools.py (eqs. (2.25)-(2.27)).
"""
from __future__ import annotations
import argparse, json, os, platform, sys, time
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from s3_chain import (SECTORS, SECTOR_TWIST, SECTOR_CLASS, E, SIGMA, R, digits, draw_couplings, hamiltonian,
                      sector_basis, block, sector_dimension_formula)
from spectral_tools import (rbar_ensemble, R_REF, unfold, surmise, amplitudes, sff_connected, cross_sff,
                            corr_matrix, sff_rmt)

HERE = os.path.dirname(os.path.abspath(__file__))
LABEL = {"e1": "([e],1)", "esgn": "([e],sgn)", "e2": "([e],2)", "t+": "([t],+)", "t-": "([t],$-$)",
         "c1": "([c],1)", "cw": r"([c],$\omega$)", "cw2": r"([c],$\omega^2$)"}
PLAIN = {"e1": "([e],1)", "esgn": "([e],sgn)", "e2": "([e],2)", "t+": "([t],+)", "t-": "([t],-)",
         "c1": "([c],1)", "cw": "([c],omega)", "cw2": "([c],omega^2)"}
# values quoted in the manuscript (Table 1)
PAPER_R = {"e1": (0.5257, 0.0010), "esgn": (0.5262, 0.0010), "e2": (0.5289, 0.0008), "t+": (0.5298, 0.0006),
           "t-": (0.5288, 0.0006), "c1": (0.5294, 0.0008), "cw": (0.5956, 0.0007), "cw2": (0.5956, 0.0007)}
L = 8                          # chain length; --L changes it, and a stored ensemble sets it
PAIR = ("cw", "cw2")
DEFAULT_SEED_BASE = 20261103


def make_ensemble(nreal, seed_base, path):
    """full diagonalization of the eight sectors for nreal realizations; realization r uses
    default_rng(seed_base + r), the same couplings for the three twists (so sectors are paired)."""
    dig = digits(L)
    bases = {s: sector_basis(L, s, dig) for s in SECTORS}
    dims = {s: bases[s].shape[1] for s in SECTORS}
    assert all(dims[s] == sector_dimension_formula(L, s) for s in SECTORS)
    spec = {s: np.empty((nreal, dims[s])) for s in SECTORS}
    coup = {"J": np.empty((nreal, L)), "f": np.empty((nreal, L)), "J2": np.empty((nreal, L - 2)), "J3": np.empty((nreal, L - 2))}
    t0 = time.time()
    for r in range(nreal):
        J, f, J2, J3 = draw_couplings(L, np.random.default_rng(seed_base + r))
        coup["J"][r], coup["f"][r], coup["J2"][r], coup["J3"][r] = J, f, J2, J3
        for h in (E, SIGMA, R):
            H = hamiltonian(L, h, J, f, J2, J3, dig)
            for s in SECTORS:
                if SECTOR_TWIST[s] == h:
                    spec[s][r] = np.linalg.eigvalsh(block(H, bases[s]))
        if r == 0 or (r + 1) % 10 == 0:
            el = time.time() - t0
            print(f"  realization {r + 1}/{nreal}  {el:.0f} s elapsed, {el / (r + 1) * nreal:.0f} s projected", flush=True)
    np.savez_compressed(path, L=L, nreal=nreal, seed_base=seed_base, sectors=np.array(SECTORS),
                        dims=np.array([dims[s] for s in SECTORS]), numpy_version=np.__version__,
                        platform=platform.platform(), **{f"spec_{s}": spec[s] for s in SECTORS},
                        **{f"coup_{k}": v for k, v in coup.items()})
    return spec, seed_base


JUNE_CODE = {"e1": "e_1", "esgn": "e_sgn", "e2": "e_2", "t+": "t_+", "t-": "t_-", "c1": "c_1", "cw": "c_w", "cw2": "c_w2"}


def load_npz(path):
    """Read a stored ensemble in either of two formats and return (spec, seed_base, description, L).
      * the format written by this script: arrays spec_<sector> of shape (nreal, dim), plus L, seed_base, ...;
      * the format of the June 2026 run (S3_rep_category_L8_spectra.npz): one array '<code>__<r>' per sector
        and realization, with codes e_1, e_sgn, e_2, t_+, t_-, c_1, c_w, c_w2 and r = 0, 1, ...; that file stores
        no seeds (its manifest gives seed = 20240601 + 1000 L + r, i.e. seed base 20248601 at L = 8)."""
    d = np.load(path, allow_pickle=False)
    if "spec_e1" in d.files:
        spec = {s: d[f"spec_{s}"] for s in SECTORS}
        seed_base = int(d["seed_base"])
        return spec, seed_base, f"seed base {seed_base}, numpy {d['numpy_version']}, {d['platform']}", int(d["L"])
    if f"{JUNE_CODE['e1']}__0" in d.files:
        nreal = 0
        while all(f"{JUNE_CODE[s]}__{nreal}" in d.files for s in SECTORS):
            nreal += 1
        spec = {s: np.array([np.sort(d[f"{JUNE_CODE[s]}__{r}"]) for r in range(nreal)]) for s in SECTORS}
        n_e = spec["e1"].shape[1] + spec["esgn"].shape[1] + 2 * spec["e2"].shape[1]
        L_file = int(round(np.log(n_e) / np.log(3)))
        if 3 ** L_file != n_e:
            sys.exit(f"ERROR: {path} does not hold complete sector spectra (dim H_e = {n_e})")
        return spec, None, f"June 2026 format, {nreal} realizations, no seeds stored", L_file
    sys.exit(f"ERROR: {path} is not an ensemble in a known format")


class Tee:
    """writes everything printed to the terminal and to a log file as well"""
    def __init__(self, *streams):
        self.streams = streams
    def write(self, text):
        for s in self.streams:
            s.write(text)
    def flush(self):
        for s in self.streams:
            s.flush()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nreal", type=int, default=120, help="number of disorder realizations, at least 2 (default 120)")
    ap.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE,
                    help="realization r uses numpy default_rng(SEED_BASE + r) (default 20261103)")
    ap.add_argument("--L", type=int, default=8, help="chain length, at least 4 (default 8; ignored with --from-npz)")
    ap.add_argument("--fraction", type=float, default=0.7, help="central fraction of each spectrum used (default 0.7)")
    ap.add_argument("--degree", type=int, default=8, help="degree of the unfolding polynomial (default 8)")
    ap.add_argument("--times", type=float, nargs=3, default=[0.02, 3.0, 80], metavar=("TMIN", "TMAX", "N"),
                    help="N log-spaced times from TMIN to TMAX, in Heisenberg times (default 0.02 3 80)")
    ap.add_argument("--band", type=float, nargs=2, default=[0.1, 1.0], metavar=("TMIN", "TMAX"),
                    help="time band over which C_AB is averaged, containing a grid time (default 0.1 1)")
    ap.add_argument("--from-npz", default=None, help="analyse a stored ensemble instead of computing one")
    ap.add_argument("--out", default=None, help="output folder (overrides the defaults described above)")
    a = ap.parse_args()
    global L
    # every value must be in range, as in the Mathematica driver, checked before any work: at L = 3 the sector
    # ([e],sgn) has too few levels, connected averages need two realizations, and C_AB needs a grid time in the band
    bad = []
    if not a.from_npz and a.nreal < 2:
        bad.append("--nreal must be at least 2")
    if not a.from_npz and a.L < 4:
        bad.append("--L must be at least 4")
    if not 0 < a.fraction <= 1:
        bad.append("--fraction must lie in (0, 1]")
    if a.degree < 1:
        bad.append("--degree must be at least 1")
    if not (np.all(np.isfinite(a.times)) and 0 < a.times[0] < a.times[1]):
        bad.append("--times needs 0 < TMIN < TMAX")
    if not (a.times[2] == int(a.times[2]) and a.times[2] >= 2):
        bad.append("--times needs a whole number N >= 2 of times")
    if not (np.all(np.isfinite(a.band)) and a.band[0] < a.band[1]):
        bad.append("--band needs TMIN < TMAX")
    if not bad:
        grid = np.logspace(np.log10(a.times[0]), np.log10(a.times[1]), int(a.times[2]))
        if not np.any((grid >= a.band[0]) & (grid <= a.band[1])):
            bad.append(f"--band {a.band[0]:g} {a.band[1]:g} contains no time of the grid --times "
                       f"{a.times[0]:g} {a.times[1]:g} {int(a.times[2])}")
    if bad:
        sys.exit("ERROR: " + ". ".join(bad))
    frac, deg, band = a.fraction, a.degree, tuple(a.band)
    tmin, tmax, ntimes = a.times[0], a.times[1], int(a.times[2])
    # output folder: results/ for the default full run, and a folder named after the non-default settings for
    # anything else, so that no other run overwrites the full-run outputs
    analysis_tags = ([f"frac{frac:g}"] if frac != 0.7 else []) + ([f"deg{deg}"] if deg != 8 else []) + \
        ([f"t{tmin:g}-{tmax:g}-{ntimes}"] if (tmin, tmax, ntimes) != (0.02, 3.0, 80) else []) + \
        ([f"band{band[0]:g}-{band[1]:g}"] if band != (0.1, 1.0) else [])
    if a.from_npz:
        if not os.path.isfile(a.from_npz):
            sys.exit(f"ERROR: ensemble file {a.from_npz} not found")
        stem = os.path.splitext(os.path.basename(a.from_npz))[0]
        res = os.path.dirname(os.path.abspath(a.from_npz))
        if os.path.samefile(res, HERE) or analysis_tags:      # never write figures into the code folder itself
            res = os.path.join(HERE, "_".join(["results_" + stem] + analysis_tags))
    else:
        L = a.L
        tags = ([f"n{a.nreal}"] if (a.nreal != 120 or a.seed_base != DEFAULT_SEED_BASE) else []) + \
            ([f"seed{a.seed_base}"] if a.seed_base != DEFAULT_SEED_BASE else []) + ([f"L{L}"] if L != 8 else []) + analysis_tags
        res = os.path.join(HERE, "_".join(["results"] + tags))
    if a.out:
        res = os.path.abspath(a.out)
    os.makedirs(res, exist_ok=True)
    log = open(os.path.join(res, "redraw_log.txt" if a.from_npz else "run_log.txt"), "w")
    sys.stdout = Tee(sys.__stdout__, log)
    t0 = time.time()
    if a.from_npz:
        spec, seed_base, description, L = load_npz(a.from_npz)
        print(f"=== stored ensemble {os.path.basename(a.from_npz)} ({description}), L = {L}")
    else:
        print(f"=== S3 qutrit chain, Section 2.4.3: L = {L}, twists e, sigma, r, {a.nreal} realizations, seed base {a.seed_base}, "
              f"numpy {np.__version__}, {platform.platform()}")
        spec, seed_base = make_ensemble(a.nreal, a.seed_base, os.path.join(res, f"ensemble_L{L}_n{a.nreal}.npz"))
    nreal = spec["e1"].shape[0]
    print(f"ensemble ready after {time.time() - t0:.0f} s")
    print(f"settings: central fraction {frac:g}, unfolding degree {deg}, {ntimes} times from {tmin:g} to {tmax:g}, "
          f"C_AB band {band[0]:g} <= t <= {band[1]:g}")
    summary = {"settings": {"L": L, "realizations": nreal, "seed_base": seed_base, "fraction": frac, "unfolding_degree": deg,
                            "corr_band": list(band), "times": [tmin, tmax, ntimes]}}

    # ---- 1. dimensions ("Each twisted Hilbert space has dimension 3^8 = 6561, and the multiplicity spaces range
    #         in dimension from 1093 to 3281.")
    dims = {s: spec[s].shape[1] for s in SECTORS}
    per_twist = {"e": dims["e1"] + dims["esgn"] + 2 * dims["e2"], "sigma": dims["t+"] + dims["t-"],
                 "r": dims["c1"] + dims["cw"] + dims["cw2"]}
    print("\n--- dimensions of the multiplicity spaces:", {PLAIN[s]: dims[s] for s in SECTORS})
    print(f"dimension of each twisted Hilbert space (([e],2) counted twice): {per_twist},  3^{L} = {3 ** L}")
    print(f"smallest and largest block: {min(dims.values())}, {max(dims.values())}   (manuscript, L = 8: 1093 to 3281)")
    summary["dimensions"] = {PLAIN[s]: dims[s] for s in SECTORS}; summary["twisted_hilbert_space_dims"] = per_twist

    # ---- 2. Table 1
    tab = {}
    print(f"\n--- Table 1: <r> from the raw gaps of the central {100 * frac:g}%, mean and standard error over realizations"
          + ("" if (L, frac) == (8, 0.7) else "  (the manuscript's values are for L = 8 and the central 70%)"))
    print(f"{'sector':14s}{'dim':>6s}{'<r>':>9s}{'SEM':>8s}{'paper':>9s}{'SEM':>8s}{'(this-paper)/sigma':>20s}{'(this-RMT)/SEM':>16s}  class")
    for s in SECTORS:
        m, sem, _, _ = rbar_ensemble(spec[s], frac)
        p, ps = PAPER_R[s]
        tab[s] = {"rbar": m, "sem": sem, "dim": dims[s], "class": SECTOR_CLASS[s], "paper": p, "paper_sem": ps}
        print(f"{PLAIN[s]:14s}{dims[s]:6d}{m:9.4f}{sem:8.4f}{p:9.4f}{ps:8.4f}{(m - p) / np.hypot(sem, ps):20.1f}"
              f"{(m - R_REF[SECTOR_CLASS[s]]) / sem:16.1f}  {SECTOR_CLASS[s]}")
    json.dump({PLAIN[s]: v for s, v in tab.items()}, open(os.path.join(res, "table1_gap_ratios.json"), "w"), indent=1)
    summary["table1"] = {PLAIN[s]: v for s, v in tab.items()}
    fig, ax = plt.subplots(figsize=(8, 4.6))
    x = np.arange(len(SECTORS))
    ax.errorbar(x, [tab[s]["rbar"] for s in SECTORS], [tab[s]["sem"] for s in SECTORS], fmt="o", color="tab:blue", label="this run")
    ax.errorbar(x + 0.15, [PAPER_R[s][0] for s in SECTORS], [PAPER_R[s][1] for s in SECTORS], fmt="s", mfc="none",
                color="tab:gray", label="Table 1 of the manuscript")
    for name, val in (("Poisson", R_REF["Poisson"]), ("GOE", R_REF["GOE"]), ("GUE", R_REF["GUE"])):
        ax.axhline(val, color="k", ls="--", lw=0.8); ax.text(-0.45, val + 0.004, name, fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels([LABEL[s] for s in SECTORS], rotation=30)
    ax.set_ylim(0.35, 0.65); ax.set_ylabel(r"$\langle r\rangle$")
    ax.set_title(f"adjacent-gap ratio, L = {L}, {nreal} realizations, central {100 * frac:g}%, raw gaps", fontsize=10.5)
    ax.legend(fontsize=9, loc="center right"); fig.tight_layout()
    fig.savefig(os.path.join(res, "table1_gap_ratios.png"), dpi=150); plt.close(fig)

    # ---- 3. Theta: Spec H_([c],omega) = Spec H_([c],omega^2), eq. (2.22)
    theta = float(np.abs(spec["cw"] - spec["cw2"]).max())
    print(f"\n--- Theta pair: max |E_n(([c],omega)) - E_n(([c],omega^2))| over all levels and realizations = {theta:.2e}")
    summary["theta_max_deviation"] = theta

    # ---- 4. Figure 3, in the style of the manuscript: one tab10 colour per sector (C0..C7 in the order of
    #         SECTORS), 50 bins on [0, 3.5], black surmises (Poisson dotted, GOE dashed, GUE dash-dotted)
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=True, sharey=True)
    sgrid = np.linspace(0, 3.5, 400)
    for k, (ax, s) in enumerate(zip(axes.flat, SECTORS)):
        gaps = np.concatenate([np.diff(unfold(E_, deg, frac)) for E_ in spec[s]])
        ax.hist(gaps, bins=50, range=(0, 3.5), density=True, alpha=0.6, color=f"C{k}", label=f"data ({len(gaps):,})")
        for c, ls in (("Poisson", ":"), ("GOE", "--"), ("GUE", "-.")):
            ax.plot(sgrid, surmise(sgrid, c), color="k", ls=ls, lw=1.0, label=c)
        ax.set_title(f"{LABEL[s]}  [{SECTOR_CLASS[s]}]   $\\langle r\\rangle$={tab[s]['rbar']:.3f}", fontsize=10.5)
        ax.legend(fontsize=7.5, loc="upper right"); ax.set_xlim(0, 3.5); ax.set_ylim(0, 1.05)
    for ax in axes[1]:
        ax.set_xlabel("s  (unfolded spacing)")
    for ax in axes[:, 0]:
        ax.set_ylabel("P(s)")
    fig.suptitle(f"Unfolded level-spacing distributions at L = {L}  ({nreal} realizations)", fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(res, "fig3_spacing.png"), dpi=150); plt.close(fig)

    # ---- 5. Figures 4 and 5, with the quantitative content of their captions
    ts = np.logspace(np.log10(tmin), np.log10(tmax), ntimes)
    Z, Ns = {}, {}
    for s in SECTORS:
        Z[s], Ns[s] = amplitudes(spec[s], ts, frac, deg)
    K = {s: sff_connected(Z[s], Ns[s]) for s in SECTORS}
    X = cross_sff(Z["cw"], Ns["cw"], Z["cw2"], Ns["cw2"])
    # Figure 4, in the style of the manuscript: tab10 colours C0..C5, GOE theory black, GUE theory red dashed
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    for k, s in enumerate([s for s in SECTORS if SECTOR_CLASS[s] == "GOE"]):
        ax.loglog(ts, K[s], "o-", color=f"C{k}", ms=3.5, lw=1.0, label=LABEL[s])
    ax.loglog(ts, sff_rmt(ts, "GOE"), "k-", lw=1.8, label="GOE theory")
    ax.loglog(ts, sff_rmt(ts, "GUE"), "r--", lw=1.2, label="GUE theory")
    ax.axvline(1.0, color="gray", ls=":", lw=1); ax.set_xlabel("t  (Heisenberg time = 1)")
    ax.set_ylabel(r"connected SFF  $K_c(t)/N$")
    ax.set_title(f"Spectral form factor of the 6 GOE-class sectors (L = {L}, {nreal} realizations)", fontsize=11)
    ax.legend(fontsize=9, loc="upper left", ncol=2); fig.tight_layout()
    fig.savefig(os.path.join(res, "fig4_sff_GOE.png"), dpi=150); plt.close(fig)
    # Figure 5, as in the manuscript
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    ax.loglog(ts, K["cw"], "o-", color="tab:red", ms=3.5, lw=1.0, label=r"auto ([c],$\omega$)")
    ax.loglog(ts, K["cw2"], "s-", color="tab:purple", ms=3.5, lw=1.0, label=r"auto ([c],$\omega^2$)")
    ax.loglog(ts, X, "D-", color="tab:green", ms=3.5, lw=1.0, label=r"cross $\omega\times\omega^2$")
    ax.loglog(ts, sff_rmt(ts, "GUE"), "-", color="red", lw=1.6, label="GUE theory")
    ax.loglog(ts, sff_rmt(ts, "GOE"), "--", color="black", lw=1.2, label="GOE theory")
    ax.axvline(1.0, color="gray", ls=":", lw=1); ax.set_xlabel("t  (Heisenberg time = 1)")
    ax.set_ylabel("connected (cross-)SFF / N")
    ax.set_title(f"SFF and cross SFF of the 2 GUE-class sectors (L = {L}, {nreal} realizations)\n"
                 r"cross $\omega\times\omega^2$ tracks the auto SFF ($\Theta$ identity)", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="upper left", ncol=2, framealpha=0.95); ax.set_xlim(ts[0] * 0.9, ts[-1] * 1.1)
    fig.tight_layout(); fig.savefig(os.path.join(res, "fig5_sff_GUE.png"), dpi=150); plt.close(fig)
    bands = [(0.02, 0.05), (0.05, 0.2), (0.2, 1.0), (1.0, 3.0)]
    ratio = {s: K[s] / sff_rmt(ts, SECTOR_CLASS[s]) for s in SECTORS}
    band_mean = lambda v, a_, b_: float(v[(ts >= a_) & (ts <= b_)].mean()) if np.any((ts >= a_) & (ts <= b_)) else float("nan")
    print("\n--- K_c / K_RMT averaged over bands of t in Heisenberg times (the manuscript states that RMT is followed for t >~ 0.2)")
    print(f"{'sector':14s}{'dim':>6s}" + "".join(f"{f'{a_}-{b_}':>11s}" for a_, b_ in bands))
    for s in SECTORS:
        print(f"{PLAIN[s]:14s}{dims[s]:6d}" + "".join(f"{band_mean(ratio[s], a_, b_):11.2f}" for a_, b_ in bands))
    sel = (ts >= 0.2) & (ts <= 3.0)
    cross_rms = float(np.sqrt(np.mean((X[sel] / K["cw"][sel] - 1) ** 2))) if sel.any() else float("nan")
    print(f"cross form factor of the Theta pair against its auto SFF, rms relative difference on 0.2 <= t <= 3: {cross_rms:.2e}")
    summary["sff_ratio_bands"] = {PLAIN[s]: {f"{a_}-{b_}": band_mean(ratio[s], a_, b_) for a_, b_ in bands} for s in SECTORS}
    summary["cross_vs_auto_rms_0.2to3"] = cross_rms

    # ---- 6. Figure 6
    keys, C, nt = corr_matrix(Z, Ns, ts, band)
    i, j = keys.index(PAIR[0]), keys.index(PAIR[1])
    off = [abs(C[p, q]) for p in range(8) for q in range(p + 1, 8) if {keys[p], keys[q]} != set(PAIR)]
    print(f"\n--- C_AB on {band[0]:g} <= t <= {band[1]:g} ({nt} grid times): Theta pair = {C[i, j]:.4f};  27 distinct pairs: "
          f"max |C| = {max(off):.3f}, rms = {np.sqrt(np.mean(np.square(off))):.3f}")
    print("          " + "".join(f"{PLAIN[k][:9]:>10s}" for k in keys))
    for p, k in enumerate(keys):
        print(f"{PLAIN[k]:14s}" + "".join(f"{C[p, q]:10.2f}" for q in range(8)))
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(C, vmin=-1, vmax=1, cmap="coolwarm")
    for p in range(8):
        for q in range(8):
            ax.text(q, p, f"{C[p, q]:.2f}", ha="center", va="center", fontsize=7.5, color="w" if abs(C[p, q]) > 0.5 else "k")
    ax.set_xticks(range(8)); ax.set_yticks(range(8))
    ax.set_xticklabels([LABEL[k] for k in keys], rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels([LABEL[k] for k in keys], fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.set_label("connected cross-SFF correlation", fontsize=9)
    ax.set_title(f"Cross-sector spectral correlations at L = {L} ({nreal} realizations)\n"
                 r"independent blocks $\to 0$;  $\Theta$ pair $(\omega,\omega^2)\to 1$", fontsize=10.5)
    fig.tight_layout(); fig.savefig(os.path.join(res, "fig6_correlation_matrix.png"), dpi=150); plt.close(fig)
    summary["C"] = {"sectors": [PLAIN[k] for k in keys], "matrix": C.tolist(), "theta_pair": float(C[i, j]),
                    "max_off_diagonal": float(max(off)), "rms_off_diagonal": float(np.sqrt(np.mean(np.square(off)))), "n_times": nt}
    json.dump(summary, open(os.path.join(res, "numerics_summary.json"), "w"), indent=1)
    print(f"\nwrote the figures and numerics_summary.json to {os.path.basename(res)}/;  total time {time.time() - t0:.0f} s")
    sys.stdout = sys.__stdout__; log.close()


if __name__ == "__main__":
    main()
