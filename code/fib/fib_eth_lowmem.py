"""
Memory-frugal, large-L ETH pipeline for the untwisted Fibonacci golden chain.

Redesign of eth_untwisted.py so that L=30 (dim ~ 1.86e6, largest momentum
block nb ~ 62k) fits the plentiful 750 GB Alliance nodes instead of needing
the scarce 6 TB tier.  Two changes vs the original pipeline:

1.  ENUMERATION  -- fib_states_recursive(L) generates only the ~phi^L
    admissible configs by DP, replacing fib_states_int's O(2^L) product().

2.  MOMENTUM BASIS kept SPARSE.  The original momentum_basis materialises a
    dense (N, nb) array (1.85 TB at L=30) and wtau_split / block_record make
    more of them.  Here the momentum basis V is a sparse (N, nb) matrix
    (~N nonzeros, 30 MB), and every block operator is formed by sparse
    products  Hb = V^H (H V)  into an nb x nb array.  Because distinct
    translation orbits have DISJOINT support, Hb (and each O_mom = V^H O V)
    is itself sparse in the orbit basis; only the eigenvectors U and the
    energy-basis observable are dense nb x nb (~60 GB each).

3.  W-SECTOR SPLIT via FLUX LABELLING, never building the full W_tau (which
    costs ~3.6 TB / 85 h at L=30: nnz(W) grows x5.83 per +2 L).  Since
    [H, W_tau] = 0, every block eigenstate psi is a pure W_tau-eigenstate, so
    its charge is  w = (W psi)_{x*} / psi_{x*}  for a SINGLE component x*.
    Choosing x* = the lowest-tau-weight config in the block (n_free ~ L/2, so
    one wtau_row of ~2^15 entries at L=30) and one gather over U gives the
    charge label of every eigenstate at once.  Validated to 1e-14 against the
    exact V^H W V split at L=16,18 (flux_label_validate.py).

Observables {eps, c2, j} are the same orbit-averaged TL words as eth_lib; the
matrix elements are computed in the energy basis as U^H (V^H O V) U.

Usage (one momentum block per task is the natural job-array unit):
    python fib_eth_lowmem.py L lam [k1 k2 ...]     # default: all k, save npz
Requires fib_chain_L22.py (numba kernels) + scinum_lowmem.py on the path.
"""

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
from math import sqrt
from numba import njit
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fib_chain_L22 import build_H_sparse, wtau_row          # tested numba kernels

phi = (1 + sqrt(5)) / 2
Fttt = np.array([[1 / phi, 1 / sqrt(phi)],
                 [1 / sqrt(phi), -1 / phi]], dtype=np.float64)

try:
    from scinum_lowmem import eigh_lowmem as _eigh_lp64   # MRRR driver, low workspace
except Exception:                                         # fallback
    def _eigh_lp64(A, overwrite_a=True, check_finite=False):
        return la.eigh(A, overwrite_a=overwrite_a, check_finite=check_finite,
                       driver="evr")


def eigh_lowmem(A, overwrite_a=True, check_finite=False):
    """Full complex-Hermitian eigendecomposition (w, V).

    Dispatches to an ILP64 zheevr (64-bit-int LAPACK) when the block dim exceeds
    the LP64 int32 index-overflow wall N > 46341 (N**2 > 2**31), which segfaults
    an LP64 LAPACK; below that, the LP64 MRRR path. Both validated to machine
    precision (see ilp64_eigh.py / ilp64_zheevr_probe.py). L<=29 -> LP64,
    L>=30 -> ILP64."""
    if A.shape[0] > 46341:
        from ilp64_eigh import eigh_ilp64
        return eigh_ilp64(A, overwrite_a=overwrite_a)
    return _eigh_lp64(A, overwrite_a=overwrite_a, check_finite=check_finite)


# ----------------------------------------------------------------- basis

def fib_states_recursive(L):
    """Sorted int64 array of admissible periodic link configs (no two adjacent
    vacua, cyclically; bit i = link i, 1 = tau).  O(D) via DP on the first bit.

    Two open-chain families by the first/last bit, then close the ring:
    enumerate all L-bit strings with no two adjacent 0s on a LINE, then keep
    those whose wrap pair (bit L-1, bit 0) is not (0,0).  The line count is
    Fib(L+2); the ring keeps Lucas(L) of them.
    """
    # DP: grow admissible OPEN strings as ints, bucketed by the last bit.
    strings_by_last = {0: [0], 1: [1]}          # length-1 strings: "0","1"
    length = 1
    while length < L:
        new0, new1 = [], []
        for s in strings_by_last[0]:            # last bit 0 -> next must be 1
            new1.append(s | (1 << length))
        for s in strings_by_last[1]:            # last bit 1 -> next 0 or 1
            new0.append(s)                      # append 0 (no new bit set)
            new1.append(s | (1 << length))
        strings_by_last = {0: new0, 1: new1}
        length += 1
    out = []
    for s in strings_by_last[0]:                # last bit 0 -> wrap ok iff bit0==1
        if s & 1:
            out.append(s)
    for s in strings_by_last[1]:                # last bit 1 -> wrap always ok
        out.append(s)
    return np.array(sorted(out), dtype=np.int64)


@njit(cache=True)
def _search(arr, val):
    lo, hi = 0, arr.size
    while lo < hi:
        mid = (lo + hi) >> 1
        if arr[mid] < val:
            lo = mid + 1
        else:
            hi = mid
    if lo < arr.size and arr[lo] == val:
        return lo
    return -1


@njit(cache=True)
def _build_ej(L, j, states, F):
    """COO of the TL generator e_j (real sparse), e_j^2 = phi e_j.
    Mirrors fib_syk_lib.Chain.e: acts at link j conditioned on links j+-1."""
    D = states.size
    jm = (j - 1) % L
    jp = (j + 1) % L
    cap = 2 * D
    rows = np.empty(cap, np.int64)
    cols = np.empty(cap, np.int64)
    vals = np.empty(cap, np.float64)
    n = 0
    w0 = F[0, 0]        # F[vac,vac] = 1/phi
    w1 = F[1, 0]        # F[tau,vac] = 1/sqrt(phi)
    for c in range(D):
        s = states[c]
        a = (s >> jm) & 1
        cc = (s >> jp) & 1
        if a != cc:
            continue
        if a == 0:                              # vacuum neighbours: diagonal phi
            rows[n] = c; cols[n] = c; vals[n] = phi; n += 1
        else:                                   # tau neighbours: phi |w><w| block
            bj = (s >> j) & 1
            wj = w0 if bj == 0 else w1
            for bnew in range(2):
                t = (s & ~(np.int64(1) << j)) | (np.int64(bnew) << j)
                ti = _search(states, t)
                if ti >= 0:
                    wb = w0 if bnew == 0 else w1
                    rows[n] = ti; cols[n] = c; vals[n] = phi * wb * wj; n += 1
    return rows[:n], cols[:n], vals[:n]


def e_gen(L, j, states):
    r, c, v = _build_ej(L, j, states, Fttt)
    D = states.size
    return sp.csr_matrix((v, (r, c)), shape=(D, D))


def untwisted_observables(L, states):
    """{eps, c2, j} as sparse complex csr on the full space (= eth_lib)."""
    D = states.size
    P = [(e_gen(L, i, states) / phi).astype(complex) for i in range(L)]
    O_eps = sp.csr_matrix((D, D), dtype=complex)
    O_c2 = sp.csr_matrix((D, D), dtype=complex)
    O_j = sp.csr_matrix((D, D), dtype=complex)
    for i in range(L):
        O_eps = O_eps + P[i]
        O_c2 = O_c2 + P[i] @ P[(i + 2) % L]
        O_j = O_j + 1j * (P[i] @ P[(i + 1) % L] - P[(i + 1) % L] @ P[i])
    return {"eps": (O_eps / L).tocsr(),
            "c2": (O_c2 / L).tocsr(),
            "j": (O_j / L).tocsr()}


# --------------------------------------------------- sparse momentum basis

def momentum_orbits(states, L):
    """List of orbits (each a list of state-row-indices under cyclic shift),
    and the index dict.  Same convention as fib_chain_L22.momentum_orbits."""
    N = states.size
    idx = {int(s): i for i, s in enumerate(states)}
    visited = np.zeros(N, dtype=bool)
    orbits = []
    mask_L = (1 << L) - 1
    for si in range(N):
        if visited[si]:
            continue
        orbit = [si]
        visited[si] = True
        s = int(states[si])
        for _ in range(L - 1):
            s = ((s >> 1) | ((s & 1) << (L - 1))) & mask_L
            sj = idx[s]
            if sj == si:
                break
            orbit.append(sj)
            visited[sj] = True
        orbits.append(orbit)
    return orbits, idx


def momentum_basis_sparse(orbits, k_idx, L, N):
    """Sparse (N, nb) momentum-k basis: column per compatible orbit, entries
    exp(-i k j)/sqrt(R) at the orbit members.  Also returns, per basis row
    (state), its (column, amplitude) so psi[y] = amp[y]*U[col[y], n] is O(1)."""
    rows, cols, vals = [], [], []
    col = 0
    for orbit in orbits:
        R = len(orbit)
        if (k_idx * R) % L != 0:
            continue
        k = 2 * np.pi * k_idx / L
        for j, si in enumerate(orbit):
            rows.append(si)
            cols.append(col)
            vals.append(np.exp(-1j * k * j) / sqrt(R))
        col += 1
    if col == 0:
        return None, None, None
    V = sp.csr_matrix((vals, (rows, cols)), shape=(N, col), dtype=complex)
    # per-row (column, amplitude) lookup (each state is in exactly one column)
    Vcsr = V.tocsr()
    row_col = np.full(N, -1, dtype=np.int64)
    row_amp = np.zeros(N, dtype=complex)
    Vcoo = V.tocoo()
    row_col[Vcoo.row] = Vcoo.col
    row_amp[Vcoo.row] = Vcoo.data
    return V.tocsc(), row_col, row_amp


# ----------------------------------------------------------- flux labels

def flux_labels(V, U, states, L, row_col, row_amp):
    """W_tau charge (+phi / -1/phi) of each block eigenstate psi_n = V @ U[:,n],
    from a SINGLE cheap wtau_row at the lowest-tau-weight supported config x*."""
    supp = np.where(row_col >= 0)[0]
    pc = np.array([bin(int(states[i])).count("1") for i in supp])
    xstar = supp[np.argmin(pc)]
    nfree = int(pc.min())
    ys, vals = wtau_row(L, int(states[xstar]), states, Fttt)
    # psi_n[y] = row_amp[y] * U[row_col[y], n]; only y in the wtau_row support
    yc = row_col[ys]
    keep = yc >= 0
    ys, vals, yc = ys[keep], vals[keep], yc[keep]
    numer = (vals * row_amp[ys]) @ U[yc, :]                  # (nb,)
    denom = row_amp[xstar] * U[row_col[xstar], :]            # (nb,)
    w = numer / denom
    plus = np.abs(w - phi) < 1e-5
    minus = np.abs(w + 1 / phi) < 1e-5
    clean = bool(np.all(plus | minus))
    return np.where(plus, +1, np.where(minus, -1, 0)), clean, nfree, w


# -------------------------------------------------- per-momentum-block solve

def block_record_energybasis(E, U, O_mom, sector_idx, wfrac=0.5,
                             max_pairs=200_000, rng=None):
    """Matrix elements of one observable inside ONE W-sector of a momentum
    block, in the energy basis:  O_E = U^H O_mom U  restricted to sector rows.
    O_mom = V^H O V (sparse nb x nb).  Returns diag + windowed off-diagonals."""
    idx = np.asarray(sector_idx)
    if idx.size < 4:
        return None
    Us = U[:, idx]                                    # (nb, ns)
    M = O_mom @ Us                                    # sparse@dense -> (nb, ns)
    diag = np.real(np.einsum('ij,ij->j', Us.conj(), M))
    imag_diag = float(np.max(np.abs(np.imag(np.einsum('ij,ij->j',
                                                      Us.conj(), M)))))
    Es = E[idx]
    order = np.argsort(Es)
    Es = Es[order]
    ns = idx.size
    lo, hi = int((1 - wfrac) / 2 * ns), int((1 + wfrac) / 2 * ns)
    win = order[lo:hi]
    Ob = Us[:, win].conj().T @ M[:, win]             # (nw, nw) energy-basis block
    Ew = E[idx][win]
    iu = np.triu_indices(len(win), k=1)
    om = np.abs(Ew[iu[1]] - Ew[iu[0]])
    vals = np.asarray(Ob)[iu]
    if len(om) > max_pairs:
        rng = rng or np.random.default_rng(0)
        keep = rng.choice(len(om), max_pairs, replace=False)
        om, vals = om[keep], vals[keep]
    return {"diag": diag[order], "imag_diag": imag_diag,
            "om": om, "o2": np.abs(vals) ** 2, "window": (lo, hi)}


def solve_block(L, lam, k_idx, states, H, obs, orbits, idx_map,
                out_prefix=None, rng=None, verbose=True):
    """One momentum block: sparse basis -> Hb_mom -> eigh -> flux labels ->
    per-sector observable matrix elements.  Returns a dict of arrays."""
    t0 = time.time()
    N = states.size
    V, row_col, row_amp = momentum_basis_sparse(orbits, k_idx, L, N)
    if V is None:
        return None
    nb = V.shape[1]
    Vh = V.conj().T.tocsr()
    Hb = (Vh @ (H @ V))                                # sparse nb x nb
    Hb = Hb.toarray()
    Hb = 0.5 * (Hb + Hb.conj().T)
    leak = None
    E, U = eigh_lowmem(Hb, overwrite_a=True)
    del Hb
    sec, clean, nfree, wvals = flux_labels(V, U, states, L, row_col, row_amp)
    data = {}
    for s, lab in ((+1, "W+"), (-1, "W-")):
        sidx = np.where(sec == s)[0]
        data[f"k{k_idx}{lab}|E"] = np.sort(E[sidx].real)
        for name, O in obs.items():
            O_mom = (Vh @ (O @ V)).tocsr()             # sparse nb x nb
            rec = block_record_energybasis(E, U, O_mom, sidx, rng=rng)
            if rec is None:
                continue
            data[f"k{k_idx}{lab}|{name}|diag"] = rec["diag"]
            data[f"k{k_idx}{lab}|{name}|om"] = rec["om"]
            data[f"k{k_idx}{lab}|{name}|o2"] = rec["o2"]
            data[f"k{k_idx}{lab}|{name}|win"] = np.array(rec["window"])
    if verbose:
        npm = int((sec == +1).sum()); nmm = int((sec == -1).sum())
        print(f"  k={k_idx:2d} nb={nb} split={npm}/{nmm} "
              f"clean={clean} n_free(x*)={nfree} ({time.time()-t0:.0f}s)",
              flush=True)
    data["_meta"] = np.array([k_idx, nb, int(clean), nfree])
    return data


def run_L(L, lam, ks=None, out_prefix="eth_lowmem", save=True, per_block=False):
    """Solve the requested momentum blocks.  per_block=True writes one npz per
    k (eth_lowmem_L{L}_lam{lam}_k{k}.npz) -- the natural job-array output, so
    concurrent array tasks never collide; otherwise one combined npz."""
    t0 = time.time()
    states = fib_states_recursive(L)
    N = states.size
    print(f"=== L={L} lam={lam}: dim={N} ===", flush=True)
    H = build_H_sparse(L, lam, states).astype(complex)
    obs = untwisted_observables(L, states)
    orbits, idx_map = momentum_orbits(states, L)
    rng = np.random.default_rng(7)
    if ks is None:
        ks = [k for k in range(1, (L + 1) // 2) if 2 * k != L]
    all_data = {}
    for k in ks:
        d = solve_block(L, lam, k, states, H, obs, orbits, idx_map, rng=rng)
        if d is None:
            continue
        payload = {kk: vv for kk, vv in d.items() if kk != "_meta"}
        if per_block and save:
            fn = f"{out_prefix}_L{L}_lam{lam}_k{k}.npz"
            np.savez_compressed(fn, **payload)
            print(f"    saved -> {fn}", flush=True)
        all_data.update(payload)
    if save and not per_block:
        fn = f"{out_prefix}_L{L}_lam{lam}.npz"
        np.savez_compressed(fn, **all_data)
        print(f"  saved -> {fn} ({time.time()-t0:.0f}s)", flush=True)
    return all_data


if __name__ == "__main__":
    # Usage: fib_eth_lowmem.py L lam [k1 k2 ...]
    #   env FIB_PER_BLOCK=1 -> one npz per k (job-array mode).  A single k may
    #   also come from SLURM_ARRAY_TASK_ID if no k args are given.
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    lam = float(sys.argv[2]) if len(sys.argv) > 2 else 1.25
    ks = [int(a) for a in sys.argv[3:]] if len(sys.argv) > 3 else None
    if ks is None and os.environ.get("SLURM_ARRAY_TASK_ID"):
        ks = [int(os.environ["SLURM_ARRAY_TASK_ID"])]
    per_block = os.environ.get("FIB_PER_BLOCK", "0") == "1"
    run_L(L, lam, ks=ks, per_block=per_block)
