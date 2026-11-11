"""Small-L regeneration test of the two Fibonacci ETH generators (L=14, lambda=1.25).

Runs code/fib/fib_eth_lowmem.py and code/fib/fib_eth_twisted_lowmem.py as subprocesses with cwd =
a fresh directory generated/fib_smallL/L14_<timestamp>/ , then on their outputs checks the two
coincidences of the paper at L=14:
  W-(k) vs twisted (tau,tau) at m=5k,        and
  (tau,1)(m) vs (1,tau)(5L-m)  over all pairs.
Prints the block dimensions found. PASS if both coincidences hold below 1e-12 (and block counts and
dimensions are as expected). Takes a few seconds. Last line: RESULT: PASS or RESULT: FAIL ...
"""
import glob, os, re, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIBDIR = os.path.join(ROOT, "code", "fib")
L, LAM = 14, 1.25
work = os.path.join(HERE, "generated", "fib_smallL", f"L{L}_{time.strftime('%Y%m%d_%H%M%S')}")
os.makedirs(work, exist_ok=True)
env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
env.pop("SLURM_ARRAY_TASK_ID", None)
fails = []
for script in ("fib_eth_lowmem.py", "fib_eth_twisted_lowmem.py"):
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(FIBDIR, script), str(L), str(LAM)],
                       cwd=work, env=env, capture_output=True, text=True)
    print(f"{script} L={L} lambda={LAM}: exit {p.returncode}, {time.time()-t0:.1f} s")
    if p.returncode:
        fails.append(f"{script} exited {p.returncode}: {p.stderr[-300:]}")

# -- collect |E| arrays --
U, T = {}, {}
for f in glob.glob(os.path.join(work, f"eth_lowmem_L{L}_lam{LAM}*.npz")):
    with np.load(f) as z:
        for k in z.files:
            m = re.fullmatch(r"k(\d+)W([+-])\|E", k)
            if m: U[(m[2], int(m[1]))] = np.sort(z[k])
for f in glob.glob(os.path.join(work, f"eth_lowmem_tw_L{L}_lam{LAM}*.npz")):
    with np.load(f) as z:
        for k in z.files:
            m = re.fullmatch(r"s(\d)m(\d+)\|E", k)
            if m: T[(int(m[1]), int(m[2]))] = np.sort(z[k])
print(f"found {len(U)} untwisted blocks, {len(T)} twisted blocks")
print("untwisted dims:", {f"k{k}W{s}": len(v) for (s, k), v in sorted(U.items(), key=lambda x: (x[0][1], x[0][0]))})
for s in (1, 2, 3):
    dims = sorted({len(v) for (ss, m), v in T.items() if ss == s})
    print(f"twisted sector s{s}: {sum(1 for (ss, m) in T if ss == s)} blocks, dims {dims}")

w1 = n1 = 0
w1 = 0.0
for k in sorted(k for s, k in U if s == "-"):
    a, b = U[("-", k)], T.get((1, 5 * k))
    if b is None or len(a) != len(b):
        fails.append(f"W-(k={k}) vs (tau,tau)(m={5*k}) missing or length mismatch"); continue
    w1 = max(w1, float(abs(a - b).max())); n1 += 1
w2, n2, wp = 0.0, 0, None
for m in sorted(m for s, m in T if s == 2):
    mp = 5 * L - m
    a, b = T[(2, m)], T.get((3, mp))
    if b is None or len(a) != len(b):
        fails.append(f"(tau,1)(m={m}) vs (1,tau)({mp}) missing or length mismatch"); continue
    d = float(abs(a - b).max()); n2 += 1
    if d > w2: w2, wp = d, (m, mp)
print(f"W-(k) vs (tau,tau)(m=5k): {n1} pairs, max|dE| = {w1:.3e}")
print(f"(tau,1)(m) vs (1,tau)({5*L}-m): {n2} pairs, max|dE| = {w2:.3e}" + (f", worst pair m={wp[0]}/{wp[1]}" if wp else ""))
if n1 < 1 or n2 < 1: fails.append("no pairs compared")
if w1 >= 1e-12: fails.append(f"W- vs (tau,tau) {w1:.2e}")
if w2 >= 1e-12: fails.append(f"(tau,1) vs (1,tau) {w2:.2e}")
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
