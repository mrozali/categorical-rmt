"""Fig. 9 (sf_fig) check: U(1)_k threshold energies at fixed spin.

Evaluates E_{s,Q}/(2 pi) = |s - (QL^2 - QR^2)/(4k)| + (QL^2 + QR^2)/(4k) - 1/6   for k=1,
(QL,QR) = (0,0) and (2,1), and compares each vertex with the TikZ coordinates of Fig. sf_fig
parsed from paper/categorical_rmt.tex (y = E/(2 pi)). Also checks the kink label (3/4, 13 pi/6).
PASS if the max deviation is < 1e-5. Nothing is written. Last line: RESULT: PASS or RESULT: FAIL ...
"""
import math, os, re

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
tex = open(os.path.join(ROOT, "paper", "categorical_rmt.tex"), encoding="utf-8").read()

def e_over_2pi(s, QL, QR, k=1):
    return abs(s - (QL**2 - QR**2) / (4 * k)) + (QL**2 + QR**2) / (4 * k) - 1 / 6

def parse_path(line_start):
    """Return the list of (x, y) vertices of the \\draw path whose line begins with line_start."""
    for ln in tex.splitlines():
        if ln.strip().startswith(line_start):
            return [(float(a), float(b)) for a, b in re.findall(r"\((-?[\d.]+),(-?[\d.]+)\)", ln)]
    raise SystemExit("RESULT: FAIL could not find TikZ line " + line_start)

paths = {
    (0, 0): parse_path(r"\draw[dashed,black,thick] (-2.4,2.233333)"),
    (2, 1): parse_path(r"\draw[cbBlue,very thick] (-2.4,4.233333)"),
}
maxdev = 0.0
for Q, pts in paths.items():
    print(f"(QL,QR)={Q}: {len(pts)} vertices parsed")
    for x, y in pts:
        f = e_over_2pi(x, *Q)
        dev = abs(f - y); maxdev = max(maxdev, dev)
        print(f"   s={x:6.3f}  paper y={y:10.6f}  formula E/2pi={f:10.6f}  |dev|={dev:.2e}")

# kink label (3/4, 13 pi/6)
m = re.search(r"\\left\(\\frac(\d)(\d),\\frac\{(\d+)\\pi\}\{(\d+)\}\\right\)", tex)
fails = []
if not m:
    fails.append("kink label not found")
else:
    s_lab = int(m[1]) / int(m[2]); E_lab = int(m[3]) / int(m[4])      # E in units of pi
    f = e_over_2pi(s_lab, 2, 1)
    dev = abs(f - E_lab / 2)
    maxdev = max(maxdev, dev)
    print(f"kink label (s,E)=({s_lab}, {m[3]}pi/{m[4]}): formula E/2pi={f:.6f}, label E/2pi={E_lab/2:.6f}, |dev|={dev:.2e}")
print(f"max deviation = {maxdev:.2e}")
if len(paths[(0, 0)]) != 3 or len(paths[(2, 1)]) != 3:
    fails.append("expected 3 vertices per path")
if not maxdev < 1e-5:
    fails.append(f"max deviation {maxdev:.2e} >= 1e-5")
print("RESULT: PASS" if not fails else "RESULT: FAIL " + "; ".join(fails))
