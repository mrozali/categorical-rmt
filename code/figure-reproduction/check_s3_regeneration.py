"""Regenerates a stored S3 ensemble from its seeds and compares it with the stored file (full rerun, ~10 min).

Usage: check_s3_regeneration.py june | nov
  june: seed base 20248601, compared with data/s3/S3_rep_category_L8_spectra.npz (Table 1, Figs. 3-4)
  nov:  seed base 20261103, compared with data/s3/ensemble_L8_n120.npz         (Figs. 5-6)
Runs code/s3/python/run_paper_numerics.py --seed-base <base> --out generated/s3_regen_<name>/, which builds the
Hamiltonian of eq. (2.23) for 120 disorder realizations of eq. (2.24) and diagonalizes all eight sectors, and then
prints the largest level-by-level difference from the stored spectra, sector by sector.
PASS if every difference is < 1e-10. Regeneration relies on numpy's random stream (numpy default_rng), so the
stored files remain the permanent record. Last line: RESULT: PASS or RESULT: FAIL ...
"""
import os, subprocess, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
GEN = os.path.join(ROOT, "code", "figure-reproduction", "generated")
CODE = os.path.join(ROOT, "code", "s3", "python")
sys.path.insert(0, CODE)
sys.dont_write_bytecode = True
from run_paper_numerics import load_npz
from s3_chain import SECTORS

RUNS = {"june": (20248601, "S3_rep_category_L8_spectra.npz"), "nov": (20261103, "ensemble_L8_n120.npz")}
if len(sys.argv) != 2 or sys.argv[1] not in RUNS:
    sys.exit(__doc__)
name = sys.argv[1]
seed_base, stored = RUNS[name]
out = os.path.join(GEN, f"s3_regen_{name}")
print(f"regenerating the {name} ensemble (seed base {seed_base}) into generated/s3_regen_{name}/ ...", flush=True)
p = subprocess.run([sys.executable, "-B", os.path.join(CODE, "run_paper_numerics.py"), "--seed-base", str(seed_base),
                    "--out", out])
if p.returncode != 0:
    sys.exit("RESULT: FAIL run_paper_numerics.py exited with an error")
new, _, _, _ = load_npz(os.path.join(out, "ensemble_L8_n120.npz"))
old, _, _, _ = load_npz(os.path.join(ROOT, "data", "s3", stored))
fails, worst = [], 0.0
print(f"\nlevel-by-level comparison with data/s3/{stored}:")
for s in SECTORS:
    if new[s].shape != old[s].shape:
        fails.append(f"{s}: shape {new[s].shape} vs {old[s].shape}")
        continue
    d = float(np.abs(np.sort(new[s], axis=1) - np.sort(old[s], axis=1)).max())
    worst = max(worst, d)
    print(f"  {s:5s} {new[s].shape[0]} realizations x {new[s].shape[1]} levels: max |dE| = {d:.2e}")
print(f"largest difference {worst:.2e}")
if worst >= 1e-10:
    fails.append(f"max |dE| = {worst:.2e}")
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
