"""
Memory-frugal, large-L ETH pipeline for the TWISTED sectors of the Fibonacci
golden chain -- the twisted analog of fib_eth_lowmem.py.

Physics is the marked-defect / twisted-translation construction of
eth_twisted.py (source of truth, validated at L<=22): the twisted translation
T_tw = g0^dag tau_cyc on the M = L+1 chain generates Z_{5L} exactly
(T_tw^L = Pi_1 + e^{4pi i/5} Pi_2 + e^{-4pi i/5} Pi_3, z_i^5 = 1), so the
canonical (sector, momentum) blocks are the Z_{5L} character projectors
P_theta = (1/5L) sum_n e^{-i n theta} T_tw^n, theta = 2 pi m / (5L), with the
sector read off from m mod 5 (0 -> sector 1, 2 -> sector 2, 3 -> sector 3;
m = 1, 4 mod 5 empty).  Per block we diagonalize H_tw(lambda) and take the
T_tw-orbit-averaged observables {eps, c2, j} from eth_lib.

WHAT IS "MEMORY-FRUGAL" HERE (mirrors fib_eth_lowmem.py where possible):

  1.  ENUMERATION -- fib_states_recursive(L+1) (imported from fib_eth_lowmem)
      generates the ~phi^{L+1} admissible configs by DP, replacing Chain's
      O(2^{L+1}) fib_states_int; cfg is bit-unpacked vectorially.  Verified
      identical to fib_states_int (validate_twisted_lowmem.py).

  2.  SPARSE OPERATORS ONLY.  T_tw, H_tw(lambda) and the three observables are
      built and kept as scipy sparse (T_tw / H_tw via sparse conjugation
      iteration, exactly as twisted_sparse_stats.build_sparse; the e_j
      generators via the numba kernel of fib_eth_lowmem, byte-identical to
      Chain.e).  W_tau is NEVER built.  No dense symmetry operator is ever
      materialized.

  3.  DENSE ONLY WHERE UNAVOIDABLE.  Because T_tw = g0^dag tau_cyc MIXES basis
      states (g0 = A + A^{-1} e_0 is not a permutation), the character-projected
      block basis Q is genuinely dense (D x nb) -- there is no sparse-orbit
      structure as in the untwisted momentum basis.  This is intrinsic to the
      eth_twisted construction (we do NOT invent a new basis).  To keep the
      peak near D x nb rather than several multiples of it we:
        * sweep the random probe batch in COLUMN CHUNKS (only one chunk of the
          T-sweep vector is live at a time; the accumulator holds D x ncols),
        * form Hb = Q^H (H Q) and each O_mom = Q^H (O Q) by CHUNKED sparse
          matvecs (never a second full D x nb dense array),
        * run ONE dense eigh per block via the IMPORTED eigh_lowmem, so the
          ILP64 zheevr path kicks in automatically at nb > 46341,
        * free every D x nb / nb x nb transient as soon as it is consumed.
      The character basis is orthonormalized by CholeskyQR2 (Gram-eigenvalue
      rank cut), the memory-lean replacement for the SVD that killed the
      L = 24 gesdd workspace, exactly as in eth_twisted.py.

Output: ONE npz per (sector, momentum) block, keyed EXACTLY as eth_twisted.py
(= the twisted analog of fib_eth_lowmem's format):
    s{i}m{m}|E                       sorted eigenvalues
    s{i}m{m}|{name}|diag             diagonal matrix elements (energy order)
    s{i}m{m}|{name}|om               |E_m - E_n| for central-window pairs
    s{i}m{m}|{name}|o2               |O_mn|^2 for those pairs
    s{i}m{m}|{name}|win              (lo, hi) window indices
for name in {eps, c2, j}.

Job-array usage (one block per task is the natural unit):
    FIB_L=30 FIB_PER_BLOCK=1 SLURM_ARRAY_TASK_ID=$k \
        python fib_eth_twisted_lowmem.py            # solves block index k
    python fib_eth_twisted_lowmem.py L lam [b1 b2 ...]   # explicit block idxs
Block index b runs over the SORTED list of non-empty m (m mod 5 in {0,2,3}).
Requires fib_eth_lowmem.py, eth_lib.py, fib_syk_lib.py on the path.
"""

import numpy as np
import scipy.sparse as sp
from math import sqrt
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fib_syk_lib import Chain, phi, fib                      # base chain + helpers
from fib_eth_lowmem import fib_states_recursive, e_gen, eigh_lowmem, Fttt
from eth_lib import twisted_observables                      # exact observables

SECTOR_OF = {0: 1, 2: 2, 3: 3}          # m mod 5 -> canonical tube sector
DEFAULT_PAD = 30                         # probe columns above the sector target


# ----------------------------------------------------------- frugal chain

class ChainLM(Chain):
    """Chain(M) with the memory-frugal, large-M enumeration:
      * states from the DP recursion (fib_states_recursive), NOT O(2^M);
      * cfg bit-unpacked vectorially;
      * e_j from the numba kernel of fib_eth_lowmem (== Chain.e byte-for-byte).
    g_j is inherited unchanged (A I + A^{-1} e_j)."""

    def __init__(self, M):
        self.M = M
        self.states = fib_states_recursive(M)
        self.D = self.states.size
        self.index = {int(s): i for i, s in enumerate(self.states)}
        self.cfg = ((self.states[:, None] >> np.arange(M)[None, :]) & 1
                    ).astype(np.int8)
        self._e = {}
        self._g = {}

    def e(self, j):
        if j not in self._e:
            self._e[j] = e_gen(self.M, j, self.states)     # real csr, == Chain.e
        return self._e[j]


def tau_cyc_sparse(ch):
    """Cyclic link shift a_j -> a_{j+1} as a sparse permutation, built
    vectorially (== twisted_sparse_stats.tau_cyc_sparse, no Python state loop):
    bit b of the image = bit (b-1) mod M of the source (left-rotate by 1)."""
    M, states = ch.M, ch.states
    mask = (1 << M) - 1
    t = ((states << 1) & mask) | (states >> (M - 1))
    rows = np.searchsorted(states, t)                       # t are valid states
    D = ch.D
    return sp.csr_matrix((np.ones(D), (rows, np.arange(D))), shape=(D, D))


def build_twisted(L, lam):
    """Sparse T_tw (complex) and H_tw(lambda) (complex) on the marked-defect
    chain M = L+1 -- identical construction to twisted_sparse_stats.build_sparse
    but on the frugally-enumerated ChainLM."""
    ch = ChainLM(L + 1)
    tau = tau_cyc_sparse(ch)
    T = (ch.g(0).getH() @ tau).tocsr()
    P1 = (ch.e(1) / phi).tocsr()
    P2 = (ch.e(2) / phi).tocsr()
    h0 = (-P1 + 0.5 * lam * (P1 @ P2 + P2 @ P1)).tocsr().astype(complex)
    H = sp.csr_matrix((ch.D, ch.D), dtype=complex)
    X = h0.copy()
    Th = T.getH().tocsr()
    for n in range(L):
        H = H + X
        X = (T @ X @ Th).tocsr()
        # prune floating-point junk outside the true (local) support
        X.data[np.abs(X.data) < 1e-13 * np.abs(X.data).max()] = 0
        X.eliminate_zeros()
    H = 0.5 * (H + H.getH())
    H = H.tocsr()
    H.data[np.abs(H.data) < 1e-14 * np.abs(H.data).max()] = 0
    H.eliminate_zeros()
    return ch, T.astype(complex), H.tocsr()


# ---------------------------------------------------- block enumeration

def block_list(L):
    """Sorted list of non-empty character momenta m (m mod 5 in {0,2,3}); the
    job-array index runs over positions in this list."""
    return [m for m in range(5 * L) if m % 5 in SECTOR_OF]


# --------------------------------------------- character basis (frugal)

# Transient-size caps (bytes) for the chunked kernels below.  Every full-size
# hidden temporary (the 51516455 OOM: chunk=None -> chunk=ncols made the sweep
# hold acc + cur + part + a matvec temp = 3-4 x the D x ncols basis) is
# replaced by a chunk bounded by these, so the per-block peak is ~1x the dense
# basis Q plus O(nb^2), not a multiple of Q.
SWEEP_BYTES = 4e9          # T-sweep column chunk (cur/part transients)
CONGR_BYTES = 2e9          # A @ Q[:, chunk] transient in congruences
GRAM_ROW_BYTES = 4e9       # row-block conj copy in the chunked Gram


def _col_chunk(D, ncols, cap_bytes):
    """Column count whose D x chunk complex block stays under cap_bytes."""
    return max(16, min(ncols, int(cap_bytes / (16 * D))))


def _gram(X, cap_bytes=GRAM_ROW_BYTES):
    """X^H X in row blocks, so X.conj() is never materialized at full size
    (np .conj() is an eager full copy -- a hidden 1x-basis transient)."""
    D, n = X.shape
    rows = max(1024, int(cap_bytes / (16 * n)))
    out = np.zeros((n, n), dtype=complex)
    for r0 in range(0, D, rows):
        blk = X[r0:r0 + rows]
        out += blk.conj().T @ blk
    return out


def _right_mult_inplace(acc, W, cap_bytes=GRAM_ROW_BYTES):
    """acc[:, :k] <- acc[:, :n] @ W for W (n x k), k <= n, row-blocked IN PLACE
    (row r of the product depends only on row r of acc, so with a row-block
    temp this is alias-safe).  Returns the leading-column view acc[:, :k];
    the D x ncols buffer is reused rather than a second D x nb array."""
    n, k = W.shape
    D = acc.shape[0]
    rows = max(1024, int(cap_bytes / (16 * n)))
    for r0 in range(0, D, rows):
        r1 = min(r0 + rows, D)
        tmp = acc[r0:r1, :n] @ W
        acc[r0:r1, :k] = tmp
    return acc[:, :k]


def character_basis(T, D, L, m, ncols, rng, chunk=None):
    """Sparse-matvec character projector P_theta applied to a random probe
    batch, orthonormalized by CholeskyQR2.  Returns the dense (D x nb) block
    basis Q with orthonormal columns (a view into the single D x ncols buffer).

    The probe batch is swept in COLUMN CHUNKS: for each chunk only its own
    T-sweep vectors (D x chunk) are live, accumulating into acc (D x ncols).
    The default chunk is SWEEP_BYTES-capped -- NOT ncols (the chunk=ncols
    default made the sweep peak 3-4 x acc and OOM-killed the L=28 run).  The
    CholeskyQR2 right-factors are applied in place in acc's buffer, so the
    whole routine holds ONE D x ncols array plus bounded transients."""
    fiveL = 5 * L
    theta = 2 * np.pi * m / fiveL
    phases = np.exp(-1j * theta * np.arange(fiveL))
    if chunk is None:
        chunk = _col_chunk(D, ncols, SWEEP_BYTES)
    acc = np.empty((D, ncols), dtype=complex)
    for c0 in range(0, ncols, chunk):
        c1 = min(c0 + chunk, ncols)
        cur = (rng.standard_normal((D, c1 - c0))
               + 1j * rng.standard_normal((D, c1 - c0)))
        part = phases[0] * cur
        for n in range(1, fiveL):
            cur = T @ cur
            part += phases[n] * cur
        acc[:, c0:c1] = part
    del cur, part
    # CholeskyQR2: rank cut on Gram eigenvalues (= squared singular values)
    Gm = _gram(acc)
    wg, Vg = np.linalg.eigh(0.5 * (Gm + Gm.conj().T))
    nb = int(np.sum(wg > 1e-14 * wg[-1]))
    Q = _right_mult_inplace(acc, Vg[:, -nb:] / np.sqrt(wg[-nb:]))
    G2 = _gram(Q)
    R = np.linalg.inv(np.linalg.cholesky(0.5 * (G2 + G2.conj().T)).conj().T)
    Q = _right_mult_inplace(acc, R)
    return Q


def sparse_congruence(Q, A, chunk=None):
    """Q^H A Q for sparse A and dense Q (D x nb), formed in column chunks so a
    full second D x nb dense array (A Q) is never materialized at once.
    Q^H itself is also never materialized (Q.conj() is an eager full copy):
    each column block is out[:, c] = Q^H Y = (Y^H Q)^H with Y = A Q[:, c]
    conjugated in place, so transients stay CONGR_BYTES-bounded.
    Returns a dense nb x nb array."""
    D, nb = Q.shape
    if chunk is None:
        chunk = _col_chunk(D, nb, CONGR_BYTES)
    out = np.empty((nb, nb), dtype=complex)
    for c0 in range(0, nb, chunk):
        c1 = min(c0 + chunk, nb)
        Y = A @ Q[:, c0:c1]                    # (D, w)
        np.conj(Y, out=Y)                      # Y -> conj(Y) in place
        out[:, c0:c1] = (Y.T @ Q).conj().T     # (Y^H Q)^H = Q^H (A Q[:, c])
    return out


def _leak_norm(H, Q, Hb, chunk=None):
    """||H Q - Q Hb||_F accumulated in column chunks (the unchunked version
    held TWO extra D x nb arrays -- the largest single transient in the old
    pipeline)."""
    D, nb = Q.shape
    if chunk is None:
        chunk = _col_chunk(D, nb, CONGR_BYTES)
    num2 = 0.0
    for c0 in range(0, nb, chunk):
        c1 = min(c0 + chunk, nb)
        Yc = H @ Q[:, c0:c1] - Q @ Hb[:, c0:c1]
        num2 += float(np.linalg.norm(Yc) ** 2)
    return sqrt(num2)


# --------------------------------------------- per-block ETH matrix elements

def block_record_from_mom(E, U, O_mom, wfrac=0.5, max_pairs=200_000, rng=None):
    """ETH summaries of one observable in one (sector, momentum) block, given
    O_mom = Q^H O Q (dense nb x nb) and the block eigenpairs (E ascending, U).
    Same window / diagonal / off-diagonal conventions as eth_lib.block_record
    (which forms the full energy-basis O); here O_E is used only on the central
    window, and the diagonal via a single O_mom @ U product."""
    nb = len(E)
    MU = O_mom @ U                                          # (nb, nb)
    d = np.einsum('ij,ij->j', U.conj(), MU)
    diag = np.real(d)
    imag_diag = float(np.max(np.abs(np.imag(d)))) if nb else 0.0
    lo, hi = int((1 - wfrac) / 2 * nb), int((1 + wfrac) / 2 * nb)
    om = np.array([])
    o2 = np.array([])
    if hi - lo >= 2:
        Uw = U[:, lo:hi]
        Ow = Uw.conj().T @ (O_mom @ Uw)                    # (nw, nw) energy basis
        Ew = E[lo:hi]
        iu = np.triu_indices(hi - lo, k=1)
        om = np.abs(Ew[iu[1]] - Ew[iu[0]])
        vals = Ow[iu]
        if len(om) > max_pairs:
            rng = rng or np.random.default_rng(0)
            keep = rng.choice(len(om), max_pairs, replace=False)
            om, vals = om[keep], vals[keep]
        o2 = np.abs(vals) ** 2
    return {"diag": diag, "imag_diag": imag_diag,
            "om": om, "o2": o2, "window": (lo, hi)}


def solve_block(L, lam, m, ch, T, H, obs, pad=DEFAULT_PAD, chunk=None,
                rng=None, verbose=True, want_matrices=False):
    """One (sector i = SECTOR_OF[m%5], momentum m) block: character basis ->
    Hb -> eigh_lowmem -> per-observable ETH summaries.  Returns a dict of
    arrays keyed s{i}m{m}|...  (plus '_meta').  If want_matrices, also returns
    the raw block eigenpairs and energy-basis observable matrices under
    '_raw' (used by the validation harness; not saved)."""
    t0 = time.time()
    i = SECTOR_OF[m % 5]
    FL, FL1 = fib(L), fib(L + 1)
    n_target = (FL1 if i == 1 else FL) // L
    ncols = n_target + pad
    if rng is None:
        rng = np.random.default_rng(1_000_000 + m)
    Q = character_basis(T, ch.D, L, m, ncols, rng, chunk=chunk)
    nb = Q.shape[1]

    Hb = sparse_congruence(Q, H)
    Hb = 0.5 * (Hb + Hb.conj().T)
    leak = float(_leak_norm(H, Q, Hb) / max(np.linalg.norm(Hb), 1e-300))
    E, U = eigh_lowmem(Hb, overwrite_a=True)
    del Hb

    lab = f"s{i}m{m}"
    data = {f"{lab}|E": np.sort(E.real)}
    raw = {"E": E, "U": U, "mats": {}} if want_matrices else None
    for name, O in obs.items():
        O_mom = sparse_congruence(Q, O)
        rec = block_record_from_mom(E, U, O_mom, rng=rng)
        data[f"{lab}|{name}|diag"] = rec["diag"]
        data[f"{lab}|{name}|om"] = rec["om"]
        data[f"{lab}|{name}|o2"] = rec["o2"]
        data[f"{lab}|{name}|win"] = np.array(rec["window"])
        if want_matrices:
            raw["mats"][name] = U.conj().T @ (O_mom @ U)   # dense nb x nb
        del O_mom
    data["_meta"] = np.array([i, m, nb, int(leak < 1e-8)])
    if verbose:
        print(f"  m={m:3d} sector {i} nb={nb} leak={leak:.1e} "
              f"({time.time()-t0:.0f}s)", flush=True)
    if want_matrices:
        return data, {"E": E, "U": U, "mats": raw["mats"], "nb": nb, "leak": leak}
    return data


def run_L(L, lam, blocks=None, out_prefix="eth_lowmem_tw", save=True,
          per_block=True, chunk=None):
    """Solve the requested block indices (positions in block_list(L)).  With
    per_block=True (default) writes one npz per block -- the job-array unit, so
    concurrent tasks never collide."""
    t0 = time.time()
    ch, T, H = build_twisted(L, lam)
    obs = twisted_observables(ch, T, L)
    ms = block_list(L)
    if blocks is None:
        blocks = list(range(len(ms)))
    print(f"=== twisted L={L} lam={lam}: dim={ch.D}, {len(ms)} blocks, "
          f"solving {len(blocks)} ===", flush=True)
    all_data = {}
    dims = {1: 0, 2: 0, 3: 0}
    for b in blocks:
        m = ms[b]
        d = solve_block(L, lam, m, ch, T, H, obs, chunk=chunk)
        i, _, nb, _ = d["_meta"]
        dims[int(i)] += int(nb)
        payload = {k: v for k, v in d.items() if k != "_meta"}
        if per_block and save:
            fn = f"{out_prefix}_L{L}_lam{lam}_b{b}.npz"
            np.savez_compressed(fn, **payload)
            print(f"    saved -> {fn}", flush=True)
        all_data.update(payload)
    if save and not per_block:
        fn = f"{out_prefix}_L{L}_lam{lam}.npz"
        np.savez_compressed(fn, **all_data)
        print(f"  saved -> {fn}", flush=True)
    if set(blocks) == set(range(len(ms))):
        print(f"  sector dims {dims} (predict "
              f"({fib(L+1)}, {fib(L)}, {fib(L)}))  total {time.time()-t0:.0f}s",
              flush=True)
    return all_data


if __name__ == "__main__":
    # Usage: fib_eth_twisted_lowmem.py [L] [lam] [b1 b2 ...]
    #   env FIB_L / FIB_PER_BLOCK / SLURM_ARRAY_TASK_ID as in fib_eth_lowmem.py.
    L = int(sys.argv[1]) if len(sys.argv) > 1 else int(os.environ.get("FIB_L", 18))
    lam = float(sys.argv[2]) if len(sys.argv) > 2 else 1.25
    blocks = [int(a) for a in sys.argv[3:]] if len(sys.argv) > 3 else None
    if blocks is None and os.environ.get("SLURM_ARRAY_TASK_ID"):
        blocks = [int(os.environ["SLURM_ARRAY_TASK_ID"])]
    per_block = os.environ.get("FIB_PER_BLOCK", "1") == "1"
    run_L(L, lam, blocks=blocks, per_block=per_block)
