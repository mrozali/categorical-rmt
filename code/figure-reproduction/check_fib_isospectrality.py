"""Isospectrality coincidences at L=26, lambda=1.25 (raw sorted eigenvalues, no unfolding).

Reads <repo root>/data/fib/fib_L26_lam1.25_spectra.npz and prints
  (i)  max |dE| between untwisted W-(k) and twisted (tau,tau) at momentum m=5k, k=1..12,
  (ii) max |dE| between twisted (tau,1) at m and (1,tau) at 130-m over all 26 pairs (worst pair named).
Expected from the reference run: 1.07e-14 and 6.48e-14 (pair m=97 / 33).
PASS if every difference is < 1e-12 (machine-precision coincidence). The bound stated in the paper
(7e-14, corrected from the 2e-14 of arXiv v1) is printed next to the measured values for reference.
Last line: RESULT: PASS or RESULT: FAIL ...
"""
import os, re
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
z = np.load(os.path.join(ROOT, "data", "fib", "fib_L26_lam1.25_spectra.npz"))
L = 26
T, U = {}, {}
for k in z.files:
    m = re.match(r"T\|s(\d)m(\d+)\|E", k)
    if m:
        T[(int(m[1]), int(m[2]))] = np.sort(z[k])
    m = re.match(r"U\|k(\d+)W([+-])\|E", k)
    if m:
        U[(m[2], int(m[1]))] = np.sort(z[k])
PAPER_BOUND = 7e-14
fails = []

ks = sorted(k for s, k in U if s == "-")
print("untwisted W- blocks k =", ks)
worst1, nblk = 0.0, 0
for k in ks:
    a, b = U[("-", k)], T.get((1, 5 * k))
    if b is None or len(a) != len(b):
        fails.append(f"W-(k={k}) vs (tau,tau)(m={5*k}): missing or length mismatch")
        continue
    worst1 = max(worst1, float(abs(a - b).max())); nblk += 1
print(f"(i)  W-(k) vs (tau,tau)(m=5k), {nblk} pairs: max|dE| = {worst1:.3e}   (paper states {PAPER_BOUND:.0e})")

worst2, wp, cnt = 0.0, None, 0
for m in sorted(m for s, m in T if s == 2):
    mp = 5 * L - m
    a, b = T[(2, m)], T.get((3, mp))
    if b is None or len(a) != len(b):
        fails.append(f"(tau,1)(m={m}) vs (1,tau)({mp}): missing or length mismatch")
        continue
    d = float(abs(a - b).max()); cnt += 1
    if d > worst2:
        worst2, wp = d, (m, mp)
print(f"(ii) (tau,1)(m) vs (1,tau)({5*L}-m), {cnt} pairs: max|dE| = {worst2:.3e} at pair m={wp[0]}/{wp[1]}   (paper states {PAPER_BOUND:.0e})")

if nblk != 12: fails.append(f"expected 12 W- pairs, got {nblk}")
if cnt != 26: fails.append(f"expected 26 (tau,1)/(1,tau) pairs, got {cnt}")
if worst1 >= 1e-12: fails.append(f"(i) {worst1:.2e} >= 1e-12")
if worst2 >= 1e-12: fails.append(f"(ii) {worst2:.2e} >= 1e-12")
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
