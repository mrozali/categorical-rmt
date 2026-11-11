"""Extract the |E| arrays of the L=26, lambda=1.25 Fibonacci banks (as read by
catrmt_fig_fib_sectors_v2.py) into ONE compact npz: data/fib/fib_L26_lam1.25_spectra.npz.
Keys: 'U|<orig key>' for untwisted (k<k>W+|E, k<k>W-|E), 'T|<orig key>' for twisted (s<s>m<m>|E).
No hygiene applied here (the plotting script applies the original hygiene).

Usage: python extract_fib_L26_spectra.py <ETH_DATA_BANK_DIR>
where the directory holds eth_lowmem_L26_lam1.25.npz and eth_lowmem_tw_L26_lam1.25_b*.npz
(outputs of fib_eth_lowmem.py / fib_eth_twisted_lowmem.py at L=26). The bank is not shipped,
the compact output is shipped in data/fib/."""
import glob, os, re, sys
import numpy as np
if len(sys.argv) != 2:
    sys.exit('usage: python extract_fib_L26_spectra.py <ETH_DATA_BANK_DIR>')
D = os.path.join(sys.argv[1], '')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'fib', 'fib_L26_lam1.25_spectra.npz')
out = {}
with np.load(D + 'eth_lowmem_L26_lam1.25.npz') as z:
    for k in z.files:
        if k.endswith('|E'):
            out['U|' + k] = z[k]
for f in sorted(glob.glob(D + 'eth_lowmem_tw_L26_lam1.25_b*.npz')):
    with np.load(f) as z:
        for k in z.files:
            if k.endswith('|E'):
                assert 'T|' + k not in out, k
                out['T|' + k] = z[k]
os.makedirs(os.path.dirname(OUT), exist_ok=True)
np.savez_compressed(OUT, **out)
print(len(out), 'arrays;', os.path.getsize(OUT) / 1e6, 'MB', os.path.abspath(OUT))
