"""
ETH matrix-element toolkit for the non-integrable Fibonacci chain,
untwisted and twisted sectors (fib_eth_note.tex).

Observables: translation-orbit averages of local TL words,
    eps = (1/L) sum_i P_i                      (identity-channel density)
    c2  = (1/L) sum_i P_i P_{i+2}              (distance-2 channel correlator)
    j   = (i/L) sum_i [P_i, P_{i+1}]           (TL current; PK-even, K-odd)
with P_i = e_i / phi.  Verified structural facts (probe run, L = 10, all
1e-16): [W_tau, P_i] = 0 for each *individual* local projector (categorical
charge is superselected for all local operators, not only H); in the twisted
sector the bare local commutes with the monodromy T^L, so the T_tw-orbit
average commutes with T_tw exactly and is block-diagonal in the canonical
(sector, momentum) blocks.

Caveat recorded in the note: off-diagonal elements of eps and of the
lambda-term density are proportional (H = -L eps + lambda L b in each
eigenbasis), so eps and c2/j are the independent off-diagonal probes.

Diagnostics per (block, observable): eigenstate expectation values (diagonal
ETH), nearest-eigenstate fluctuations, off-diagonal pairs (omega, |O_mn|^2)
in a central energy window, the GOE/GUE variance ratio
Sigma^2_d / Sigma^2_od (2 for GOE, 1 for GUE), and the basis-phase-invariant
moduli ratio Gamma = <|O|^2> / <|O|>^2 (pi/2 for real/GOE-class Gaussian,
4/pi for complex/GUE-class Gaussian).
"""

import numpy as np
import scipy.sparse as sp

phi = (1 + np.sqrt(5)) / 2


# ------------------------------------------------------------- observables

def untwisted_observables(ch):
    """Zero-momentum orbit averages on the periodic golden chain
    (ch = fib_syk_lib.Chain(L), M = L links).  Returns {name: csr complex}."""
    L = ch.M
    P = [(ch.e(i) / phi).tocsr().astype(complex) for i in range(L)]
    D = ch.D
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


def twisted_observables(ch, T, L):
    """T_tw-orbit averages of the same local words on the marked-defect
    chain (M = L+1 links, L physical anyons).  Base densities sit at bulk
    links 1..3, away from the seam; the orbit supplies the consistently
    braid-dressed seam terms, exactly as in H_tw."""
    P1 = (ch.e(1) / phi).tocsr().astype(complex)
    P2 = (ch.e(2) / phi).tocsr().astype(complex)
    P3 = (ch.e(3) / phi).tocsr().astype(complex)
    seeds = {"eps": P1,
             "c2": P1 @ P3,
             "j": 1j * (P1 @ P2 - P2 @ P1)}
    Th = T.getH().tocsr()
    out = {}
    for name, h0 in seeds.items():
        X = h0.tocsr()
        acc = sp.csr_matrix((ch.D, ch.D), dtype=complex)
        for n in range(L):
            acc = acc + X
            if n < L - 1:
                X = (T @ X @ Th).tocsr()
                # prune float junk outside the true local support
                X.data[np.abs(X.data) < 1e-13 * np.abs(X.data).max()] = 0
                X.eliminate_zeros()
        out[name] = (acc / L).tocsr()
    return out


# ------------------------------------------------------- per-block records

def block_record(E, U, Vb, O_sparse, wfrac=0.5, max_pairs=200_000,
                 rng=None):
    """Matrix elements of one observable in one symmetry block.

    E, U   : eigenvalues / eigenvectors of the block Hamiltonian (ascending);
    Vb     : (D x nb) block basis (columns orthonormal);
    O_sparse: the observable on the full Hilbert space.

    Returns dict with the diagonal, and (omega, |O|^2, Re, Im) for all
    off-diagonal pairs whose mean energy lies in the central `wfrac`
    fraction of the block spectrum (subsampled to max_pairs).
    Also returns the block leakage of O (exact-zero gate)."""
    OV = O_sparse @ Vb
    leak = np.linalg.norm(OV - Vb @ (Vb.conj().T @ OV)) / \
        max(np.linalg.norm(OV), 1e-300)
    Ob = U.conj().T @ (Vb.conj().T @ OV) @ U
    nb = len(E)
    diag = np.real(np.diag(Ob))
    imag_diag = np.max(np.abs(np.imag(np.diag(Ob)))) if nb else 0.0
    lo, hi = int((1 - wfrac) / 2 * nb), int((1 + wfrac) / 2 * nb)
    idx = np.arange(lo, hi)
    om_list, o2_list = [], []
    if len(idx) >= 2:
        sub = Ob[np.ix_(idx, idx)]
        Ei = E[idx]
        iu = np.triu_indices(len(idx), k=1)
        om = np.abs(Ei[iu[1]] - Ei[iu[0]])
        vals = sub[iu]
        if len(om) > max_pairs:
            rng = rng or np.random.default_rng(0)
            keep = rng.choice(len(om), max_pairs, replace=False)
            om, vals = om[keep], vals[keep]
        om_list, o2_list = om, np.abs(vals) ** 2
    return {"diag": diag, "imag_diag": imag_diag, "leak": leak,
            "om": np.asarray(om_list), "o2": np.asarray(o2_list),
            "window": (lo, hi)}


# ------------------------------------------------------------- diagnostics

def sigma2_diag(E, diag, lo, hi):
    """Detrended diagonal variance from nearest-eigenstate differences:
    Sigma^2_d = <|O_{n+1,n+1} - O_nn|^2> / 2 over the window (kills the
    smooth O(E) trend at leading order)."""
    d = np.diff(diag[lo:hi])
    return 0.5 * np.mean(d ** 2) if len(d) else np.nan


def small_omega_stats(om, o2, n_spacings, span, c=6.0):
    """Off-diagonal variance and moduli ratio Gamma at small omega:
    pairs with omega < c * mean level spacing (span / n_spacings)."""
    if len(om) == 0 or n_spacings == 0:
        return np.nan, np.nan, 0
    cut = c * span / n_spacings
    sel = om < cut
    if sel.sum() < 10:
        return np.nan, np.nan, int(sel.sum())
    s2 = np.mean(o2[sel])
    gam = s2 / np.mean(np.sqrt(o2[sel])) ** 2
    return s2, gam, int(sel.sum())


def gamma_of_omega(om, o2, nbins=24):
    """Gamma(omega) = <|O|^2>/<|O|>^2 in log-spaced omega bins.
    pi/2 = real Gaussian (GOE class), 4/pi = complex Gaussian (GUE)."""
    pos = om > 0
    om, o2 = om[pos], o2[pos]
    if len(om) < 100:
        return np.array([]), np.array([]), np.array([])
    edges = np.geomspace(om.min() * 0.999, om.max() * 1.001, nbins + 1)
    which = np.digitize(om, edges) - 1
    centers, gam, cnt = [], [], []
    for b in range(nbins):
        sel = which == b
        if sel.sum() < 30:
            continue
        centers.append(np.sqrt(edges[b] * edges[b + 1]))
        gam.append(np.mean(o2[sel]) / np.mean(np.sqrt(o2[sel])) ** 2)
        cnt.append(sel.sum())
    return np.array(centers), np.array(gam), np.array(cnt)


def f2_profile(om, o2, rho, nbins=40):
    """Spectral-function profile |f(Ebar, omega)|^2 ~ rho * <|O_mn|^2>
    binned in omega (rho = local density of states of the block window)."""
    pos = om > 0
    om, o2 = om[pos], o2[pos]
    if len(om) < 100:
        return np.array([]), np.array([])
    edges = np.geomspace(max(om.min(), 1e-12), om.max() * 1.001, nbins + 1)
    which = np.digitize(om, edges) - 1
    centers, prof = [], []
    for b in range(nbins):
        sel = which == b
        if sel.sum() < 20:
            continue
        centers.append(np.sqrt(edges[b] * edges[b + 1]))
        prof.append(rho * np.mean(o2[sel]))
    return np.array(centers), np.array(prof)


def goe_calibration(n, seeds=24, obs="proj"):
    """Finite-size calibration of (ratio, Gamma) on synthetic GOE at matched
    block size n.  Measured (quoted in the note): n=143: 2.11(10)/1.559(9);
    n=232: 2.07(6)/1.578(10); n=547: 2.04(4)/1.571(5); n=1270:
    2.01(3)/1.568(2) -- i.e. estimator bias is small at production sizes;
    the physical excess over 2 is genuine finite-size ETH correction."""
    rng2 = np.random.default_rng(999)
    rows = []
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        Hg = rng.standard_normal((n, n))
        Hg = (Hg + Hg.T) / np.sqrt(2 * n)
        E, U = np.linalg.eigh(Hg)
        if obs == "proj":
            O = np.zeros((n, n))
            O[np.arange(n // 2), np.arange(n // 2)] = 1.0
        else:
            O = rng2.standard_normal((n, n))
            O = (O + O.T) / 2
        Ob = U.T @ O @ U
        lo, hi = n // 4, 3 * n // 4
        d = np.diag(Ob)[lo:hi]
        var_nn = 0.5 * np.mean(np.diff(d) ** 2)
        idx = np.arange(lo, hi)
        iu = np.triu_indices(len(idx), k=1)
        om = np.abs(E[idx][iu[1]] - E[idx][iu[0]])
        o2 = np.abs(Ob[np.ix_(idx, idx)][iu]) ** 2
        sel = om < 6 * (E[hi] - E[lo]) / (hi - lo)
        rows.append((var_nn / np.mean(o2[sel]),
                     np.mean(o2[sel]) / np.mean(np.sqrt(o2[sel])) ** 2))
    r, g = np.array(rows).T
    return (r.mean(), r.std() / np.sqrt(seeds)), \
           (g.mean(), g.std() / np.sqrt(seeds))


def goe_pipeline_sanity(n=800, seed=1):
    """Tier-1 gate: GOE H + fixed projector observable must give
    Sigma^2_d/Sigma^2_od = 2 and Gamma = pi/2 (no physics)."""
    rng = np.random.default_rng(seed)
    Hg = rng.standard_normal((n, n))
    Hg = (Hg + Hg.T) / np.sqrt(2 * n)
    E, U = np.linalg.eigh(Hg)
    O = np.zeros((n, n))
    O[np.arange(n // 2), np.arange(n // 2)] = 1.0
    Ob = U.T @ O @ U
    lo, hi = n // 4, 3 * n // 4
    s2d = sigma2_diag(E, np.diag(Ob), lo, hi)
    idx = np.arange(lo, hi)
    iu = np.triu_indices(len(idx), k=1)
    om = np.abs(E[idx][iu[1]] - E[idx][iu[0]])
    o2 = np.abs(Ob[np.ix_(idx, idx)][iu]) ** 2
    s2od, gam, cnt = small_omega_stats(om, o2, hi - lo, E[hi] - E[lo])
    return s2d / s2od, gam
