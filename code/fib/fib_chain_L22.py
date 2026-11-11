"""
L=22 Fibonacci chain with full Z(Fib) tube-algebra resolution.
W_tau built via Numba JIT (rows in parallel via Python-level outer loop).
H built sparsely to avoid the dense-N^2 memory wall.
"""

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from numba import njit
from itertools import product
from math import sqrt
import os, time

phi = (1 + sqrt(5)) / 2
Fttt_t = np.array([[1/phi,        1/sqrt(phi)],
                   [1/sqrt(phi), -1/phi      ]], dtype=np.float64)


# --- state enumeration ---

def fib_states_int(L):
    out = []
    for bits in product((0, 1), repeat=L):
        ok = True
        for i in range(L):
            if bits[i] == 0 and bits[(i + 1) % L] == 0:
                ok = False; break
        if ok:
            v = 0
            for i, b in enumerate(bits):
                v |= b << i
            out.append(v)
    out.sort()
    return np.array(out, dtype=np.int64)


@njit(cache=True, inline='always')
def _bit(v, i):
    return (v >> i) & 1


@njit(cache=True)
def _search(arr, val):
    lo = 0; hi = arr.size
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
def wtau_row(L, x_int, states_int, F):
    """Nonzero entries in column x of W_tau: returns (row indices, values)."""
    x = np.empty(L, dtype=np.int8)
    for i in range(L):
        x[i] = _bit(x_int, i)
    free = np.empty(L, dtype=np.int64)
    n_free = 0
    for i in range(L):
        if x[i] == 1:
            free[n_free] = i; n_free += 1
    cap = 1 << n_free
    out_idx = np.empty(cap, dtype=np.int64)
    out_val = np.empty(cap, dtype=np.float64)
    count = 0
    y = np.empty(L, dtype=np.int8)
    for mask in range(cap):
        for i in range(L):
            if x[i] == 0:
                y[i] = 1
        for j in range(n_free):
            i = free[j]
            y[i] = np.int8((mask >> j) & 1)
        ok = True
        for i in range(L):
            ip1 = i + 1
            if ip1 == L: ip1 = 0
            if y[i] == 0 and y[ip1] == 0:
                ok = False; break
        if not ok: continue
        y_int = 0
        for i in range(L):
            y_int |= np.int64(y[i]) << i
        sb_i = _search(states_int, y_int)
        if sb_i < 0: continue
        val = 1.0
        for i in range(L):
            ip = i - 1
            if ip < 0: ip = L - 1
            ap = x[ip]; ac = x[i]; bp = y[ip]; bc = y[i]
            if ap == 1 and bc == 1:
                val *= F[ac, bp]
        out_idx[count] = sb_i
        out_val[count] = val
        count += 1
    return out_idx[:count], out_val[:count]


def build_Wtau(L, states_int, verbose=True):
    N = states_int.size
    # warm up JIT
    _ = wtau_row(L, states_int[0], states_int, Fttt_t)
    t0 = time.time()
    rows_chunks, cols_chunks, vals_chunks = [], [], []
    nnz_total = 0
    for sa_i in range(N):
        idxs, vals = wtau_row(L, states_int[sa_i], states_int, Fttt_t)
        if len(idxs):
            rows_chunks.append(idxs)
            cols_chunks.append(np.full(len(idxs), sa_i, dtype=np.int64))
            vals_chunks.append(vals)
            nnz_total += len(idxs)
        if verbose and sa_i % 5000 == 0 and sa_i > 0:
            print(f"    {sa_i}/{N}  ({time.time()-t0:.1f}s, {nnz_total} nnz)")
    if verbose:
        print(f"  inner loop done in {time.time()-t0:.1f}s, total nnz = {nnz_total}")
    t0 = time.time()
    rows = np.concatenate(rows_chunks)
    cols = np.concatenate(cols_chunks)
    vals = np.concatenate(vals_chunks)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))
    if verbose:
        print(f"  CSR conversion: {time.time()-t0:.1f}s, density = {W.nnz / N**2:.4f}")
    return W


# --- H sparse ---

@njit(cache=True)
def _H_entries_for_state(L, x_int, states_int, F, lam, out_rows, out_cols, out_vals, count):
    # Apply -sum_i P_i^{(1)}  and  +lam sum_i P_i^{(1)} P_{i+1}^{(1)}
    x = np.empty(L, dtype=np.int8)
    for i in range(L):
        x[i] = _bit(x_int, i)
    si = _search(states_int, x_int)
    # 1-site terms
    for i in range(L):
        ip = i - 1
        if ip < 0: ip = L - 1
        ip1 = i + 1
        if ip1 == L: ip1 = 0
        ap = x[ip]; an = x[ip1]; ac = x[i]
        if ap == 0 and an == 0:
            out_rows[count] = si; out_cols[count] = si
            out_vals[count] = -1.0
            count += 1
        elif ap == 0 or an == 0:
            pass
        else:
            for ac_new in range(2):
                y_int = (x_int & ~(np.int64(1) << i)) | (np.int64(ac_new) << i)
                sj = _search(states_int, y_int)
                if sj >= 0:
                    out_rows[count] = sj; out_cols[count] = si
                    out_vals[count] = -F[ac_new, 0] * F[ac, 0]
                    count += 1
    # NNN term: P_i P_{i+1}, apply P_i then P_{i+1}
    if lam != 0.0:
        # temp store: up to 2 intermediate states with amplitudes
        i_states = np.empty(2, dtype=np.int64)
        i_amps = np.empty(2, dtype=np.float64)
        for i in range(L):
            ip = i - 1
            if ip < 0: ip = L - 1
            ip1 = i + 1
            if ip1 == L: ip1 = 0
            ap = x[ip]; an = x[ip1]; ac = x[i]
            n_int = 0
            if ap == 0 and an == 0:
                i_states[0] = x_int; i_amps[0] = 1.0
                n_int = 1
            elif ap == 0 or an == 0:
                continue
            else:
                for ac_new in range(2):
                    y_int = (x_int & ~(np.int64(1) << i)) | (np.int64(ac_new) << i)
                    i_states[n_int] = y_int
                    i_amps[n_int] = F[ac_new, 0] * F[ac, 0]
                    n_int += 1
            # apply P_{i+1}^(1) to each intermediate
            im1 = ip1
            for k in range(n_int):
                y_int = i_states[k]; amp = i_amps[k]
                yp = _bit(y_int, im1 - 1 if im1 > 0 else L - 1)
                yn = _bit(y_int, im1 + 1 if im1 < L - 1 else 0)
                yc = _bit(y_int, im1)
                if yp == 0 and yn == 0:
                    sj = _search(states_int, y_int)
                    if sj >= 0:
                        out_rows[count] = sj; out_cols[count] = si
                        out_vals[count] = lam * amp
                        count += 1
                elif yp == 0 or yn == 0:
                    pass
                else:
                    for yc_new in range(2):
                        z_int = (y_int & ~(np.int64(1) << im1)) | (np.int64(yc_new) << im1)
                        sj = _search(states_int, z_int)
                        if sj >= 0:
                            out_rows[count] = sj; out_cols[count] = si
                            out_vals[count] = lam * amp * F[yc_new, 0] * F[yc, 0]
                            count += 1
    return count


def build_H_sparse(L, lam, states_int):
    N = states_int.size
    # Upper bound on nnz per state: ~ 6L (very generous)
    cap = N * 8 * L
    rows = np.empty(cap, dtype=np.int64)
    cols = np.empty(cap, dtype=np.int64)
    vals = np.empty(cap, dtype=np.float64)
    count = 0
    # warm up JIT
    _ = _H_entries_for_state(L, states_int[0], states_int, Fttt_t, lam,
                             rows, cols, vals, 0)
    count = 0
    t0 = time.time()
    for si in range(N):
        count = _H_entries_for_state(L, states_int[si], states_int, Fttt_t, lam,
                                      rows, cols, vals, count)
    print(f"  H sparse loop: {time.time()-t0:.1f}s, raw nnz = {count}")
    H = sp.csr_matrix((vals[:count], (rows[:count], cols[:count])), shape=(N, N))
    H = 0.5 * (H + H.T)
    return H


# --- symmetry resolution ---

def momentum_orbits(states_int, L):
    N = states_int.size
    idx = {int(s): i for i, s in enumerate(states_int)}
    visited = np.zeros(N, dtype=bool)
    orbits = []
    mask_L = (1 << L) - 1
    for si in range(N):
        if visited[si]: continue
        orbit = [si]; visited[si] = True
        s = int(states_int[si])
        for _ in range(L - 1):
            # cyclic left shift by 1: take bit 0 to position L-1
            s = ((s >> 1) | ((s & 1) << (L - 1))) & mask_L
            sj = idx[s]
            if sj == si: break
            orbit.append(sj); visited[sj] = True
        orbits.append(orbit)
    return orbits, idx


def momentum_basis(orbits, k_idx, L, N):
    cols = []
    for orbit in orbits:
        R = len(orbit)
        if (k_idx * R) % L != 0: continue
        k = 2 * np.pi * k_idx / L
        v = np.zeros(N, dtype=complex)
        for j, si in enumerate(orbit):
            v[si] = np.exp(-1j * k * j) / sqrt(R)
        cols.append(v)
    return np.column_stack(cols) if cols else None


def parity_split(V, states_int, idx, L):
    def reflect(s):
        r = 0
        for i in range(L):
            r |= ((s >> i) & 1) << (L - 1 - i)
        return r
    d = V.shape[1]; N = V.shape[0]
    Pmat = np.zeros((d, d), dtype=complex)
    # precompute reflection map
    refl_map = np.empty(N, dtype=np.int64)
    for si in range(N):
        refl_map[si] = idx[reflect(int(states_int[si]))]
    for a in range(d):
        v = V[:, a]
        v_ref = np.zeros(N, dtype=complex)
        nz = np.nonzero(np.abs(v) > 1e-14)[0]
        for si in nz:
            v_ref[refl_map[si]] += v[si]
        Pmat[:, a] = V.conj().T @ v_ref
    Pmat = 0.5 * (Pmat + Pmat.conj().T)
    w, U = la.eigh(Pmat)
    return V @ U[:, w > 0], V @ U[:, w <= 0]


def wtau_split(V, W):
    Wb = V.conj().T @ (W @ V)
    Wb = 0.5 * (Wb + Wb.conj().T)
    w, U = la.eigh(Wb)
    return V @ U[:, w > 0], V @ U[:, w <= 0]


def spacings(e, drop=0.05):
    e = np.sort(np.real(e)); n = len(e); lo = int(drop * n); hi = n - lo
    if hi - lo < 5: return np.array([])
    s = np.diff(e[lo:hi]); s = s[s > 1e-9]
    return s / np.mean(s) if len(s) else s


def ratios(e, drop=0.05):
    e = np.sort(np.real(e)); n = len(e); lo = int(drop * n); hi = n - lo
    if hi - lo < 6: return np.array([])
    s = np.diff(e[lo:hi])
    if len(s) < 2: return np.array([])
    r = np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])
    return r[np.isfinite(r)]


def main(L=22, lam=1.25):
    print(f"=== Fibonacci anyon chain, L = {L}, lambda = {lam} ===")
    print("Enumerating states...")
    t0 = time.time()
    states_int = fib_states_int(L)
    N = states_int.size
    print(f"  dim = {N}  (Lucas({L})),  {time.time()-t0:.1f}s")

    print("Building H sparsely (Numba)...")
    t0 = time.time()
    H = build_H_sparse(L, lam, states_int)
    print(f"  H nnz = {H.nnz},  built in {time.time()-t0:.1f}s")

    Wpath = f"Wtau_L{L}.npz"
    if os.path.exists(Wpath):
        print(f"Loading W_tau from {Wpath} ...")
        W = sp.load_npz(Wpath)
    else:
        print("Building W_tau (Numba)...")
        W = build_Wtau(L, states_int)
        sp.save_npz(Wpath, W)
        print(f"  cached to {Wpath}")

    print("\nVerifying W_tau:")
    # Quick spot check: commutator on a few random vectors
    rng = np.random.default_rng(0)
    v = rng.standard_normal(N)
    HW = H @ (W @ v)
    WH = W @ (H @ v)
    rel = np.linalg.norm(HW - WH) / np.linalg.norm(HW + WH + 1e-30)
    print(f"  ||[H,W]v|| / ||{{H,W}}v|| (random v) = {rel:.2e}")
    # Symmetric?
    Wd = W - W.T
    sym_rel = sp.linalg.norm(Wd) / sp.linalg.norm(W)
    print(f"  symmetric: ||W - W^T|| / ||W|| = {sym_rel:.2e}")

    print("\nResolving (k, parity, W_tau) blocks ...")
    orbits, idx = momentum_orbits(states_int, L)
    results = []
    t0 = time.time()
    for k_idx in range(L):
        V = momentum_basis(orbits, k_idx, L, N)
        if V is None: continue
        if k_idx == 0 or (L % 2 == 0 and k_idx == L // 2):
            Vp, Vm = parity_split(V, states_int, idx, L)
            sub_blocks = [(Vp, +1), (Vm, -1)]
        else:
            sub_blocks = [(V, 0)]
        for Vkp, par in sub_blocks:
            if Vkp.shape[1] == 0: continue
            Vplus, Vminus = wtau_split(Vkp, W)
            for Vsub, wlab in ((Vplus, +1), (Vminus, -1)):
                if Vsub.shape[1] == 0: continue
                Hb = Vsub.conj().T @ (H @ Vsub)
                Hb = 0.5 * (Hb + Hb.conj().T)
                ev = la.eigvalsh(Hb)
                results.append(((k_idx, par, wlab), ev))
        print(f"  k={k_idx} done  ({time.time()-t0:.1f}s, {len(results)} blocks)")

    total = sum(len(e) for _, e in results)
    print(f"\n  total blocks: {len(results)}, dim sum: {total}  (N = {N})")

    s_plus = np.concatenate([spacings(e) for lab, e in results if lab[2] == +1])
    s_minus = np.concatenate([spacings(e) for lab, e in results if lab[2] == -1])
    r_plus = np.concatenate([ratios(e) for lab, e in results if lab[2] == +1])
    r_minus = np.concatenate([ratios(e) for lab, e in results if lab[2] == -1])
    print(f"  vacuum-flux: {len(s_plus)} spacings,  <r>={np.mean(r_plus):.4f}")
    print(f"  tau-flux:    {len(s_minus)} spacings, <r>={np.mean(r_minus):.4f}")
    print(f"  (GOE expects <r> = 0.5307)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    s_grid = np.linspace(0.001, 4, 500)
    GOE = (np.pi / 2) * s_grid * np.exp(-np.pi / 4 * s_grid ** 2)
    Pois = np.exp(-s_grid)
    for ax, sdata, rdata, color, title, flux in [
        (axes[0], s_plus, r_plus, "seagreen",
         "$W_\\tau = +\\varphi$  (vacuum flux)", "vacuum flux"),
        (axes[1], s_minus, r_minus, "darkorange",
         "$W_\\tau = -1/\\varphi$  ($\\tau$ flux)", "$\\tau$ flux"),
    ]:
        ax.hist(sdata, bins=60, range=(0, 4), density=True, alpha=0.55, color=color,
                label=f"{flux}  ({len(sdata)} spacings, $\\langle r\\rangle$={np.mean(rdata):.4f})")
        ax.plot(s_grid, GOE, 'b--', lw=2, label="single GOE  ($\\langle r\\rangle$=0.5307)")
        ax.plot(s_grid, Pois, 'r--', lw=1.4, label="Poisson")
        ax.set_xlim(0, 4); ax.set_ylim(0, 1.0)
        ax.set_xlabel("s"); ax.set_ylabel("P(s)")
        ax.set_title(title, fontsize=13)
        ax.legend(fontsize=10, loc="upper right")
    plt.suptitle(f"L=22 Fibonacci chain (dim={N}), $\\lambda$=1.25:  "
                 f"$P(s)$ within $(k, \\mathcal{{P}}, W_\\tau)$ blocks (5% band-edge drop)",
                 fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"fib_single_goe_L{L}.png", dpi=140, bbox_inches="tight")
    print(f"\nWrote fib_single_goe_L{L}.png")


if __name__ == "__main__":
    main(L=22)
