"""
Library for the anyonic SYK model on the Fibonacci golden chain
(anyonic_syk_notes.tex) and its tau-twisted sector.

Conventions
-----------
* Periodic chain of M tau anyons <-> M links. Link j in {0=vacuum, 1=tau},
  bit-encoded exactly as in fib_chain/fib_chain_L22.py (bit i = link i);
  admissibility: no two adjacent vacuum links.
* Anyon j (j = 1..M) sits between link (j-1) and link (j mod M).
  Link 0 is the SEAM: interactions are restricted to subsets
  S subset {1..N_phys} whose arc [min S, max S] never wraps through it
  (fundamental-domain convention; the seam link enters at most as a
  read-only boundary label).
* Untwisted model: M = N (all anyons physical).
  Twisted sector: M = N+1, the marked-defect chain -- anyon "0" between
  links M-1 and 0 is the tau defect; physical anyons are 1..N and the
  interaction subsets are literally the same as untwisted.  Since the
  defect line is topological and no interaction arc crosses the seam,
  this Hilbert-space substitution (with unchanged Hamiltonian terms)
  IS the categorical twist.  (No seam term arises, in contrast with the
  translation-invariant golden chain -- disorder breaks translation, so
  the fundamental domain is a legitimate gauge choice.)
* F-symbol F^{tau tau tau}_tau in the (vacuum, tau) = (index 0, index 1)
  basis, as in the notes (eq. after "Non-triviality..."):
      F = [[1/phi, 1/sqrt(phi)], [1/sqrt(phi), -1/phi]].
  All F-symbols with any external label = vacuum equal 1 when admissible.
* Temperley-Lieb e_j acts at link j conditioned on links j+-1 (same
  construction verified in fib_chain/exploration/golden_full.py:
  e^2 = phi e, e e' e = e for neighbours, distant commute).
  Pair identity-channel projector of anyons (j, j+1):  P^(1) = e_j / phi.
* Braid g_j = A + A^{-1} e_j, A = exp(3 i pi / 5) (Kauffman, unitary;
  braid relations verified in golden_full.py / rsos_phase2_l4.py).
  g_j exchanges anyons (j, j+1) with the "over" chirality; g_j^dagger
  is the "under" exchange.
* P^(1)_S for S = {s_1 < ... < s_k} (notes eq. (1)/(4), P_S = F_S^dag P_root F_S):
    1. braid-transport the anyons of S into the contiguous block
       [s_k - k + 1, s_k]  (compress right; each S-anyon crosses
       spectators with a fixed chirality),
    2. project the block onto collective channel vacuum via the F-move
       "zip" (an exact isometry to the fused-tree basis),
    3. transport back:  P^(1)_S = G_S^dagger P^(1)_block G_S.

  CHIRALITY (correcting Remark 1 of the notes): for non-contiguous S the
  channel projector is NOT routing-independent.  Routings fall into two
  classes by the side on which S-strands pass the spectators; the two
  resulting projectors are complex conjugates of each other,
      P^-_S = (P^+_S)^*   (both Hermitian idempotents, same rank),
  and differ at O(1).  The would-be path-independence argument fails
  because transport moves single tau strands (not the fused channel-1
  cluster) past spectators, and single-strand monodromy is nontrivial in
  both fusion channels.  Within a fixed chirality class the operator IS
  canonical (verified: the compress-right and compress-left words agree
  once their crossing sides are matched).  The model therefore comes in
  two time-reversal-conjugate chiral versions; we adopt one ('+') as the
  definition, so H(J) is complex Hermitian and TRS is broken by each
  interaction term individually (not by the disorder).

Everything is dense/sparse numpy at Lucas(M) <= a few thousand.
"""

import numpy as np
import scipy.sparse as sp
from math import sqrt
from numba import njit

phi = (1 + sqrt(5)) / 2
A_KAUFF = np.exp(3j * np.pi / 5)          # Kauffman A with d = -A^2 - A^-2 = phi
Fttt = np.array([[1 / phi, 1 / sqrt(phi)],
                 [1 / sqrt(phi), -1 / phi]], dtype=np.float64)

# Ferraggata-van Rees lasso eigenvalues for L^{tau tau}_tau
# (fib_chain/exploration/test_complex_structure.py conventions)
LAM1 = sqrt(phi) / (1 + phi)                                   # sector (1)
LAM2 = sqrt(phi) / (1 + phi * np.exp(4j * np.pi / 5))          # sector (2)
LAM3 = np.conj(LAM2)                                           # sector (3)


# ---------------------------------------------------------------- basis

def fib_states_int(M):
    """Sorted int64 array of admissible periodic link configurations.
    (= fib_chain_L22.fib_states_int; bit i = link i, 1 = tau.)"""
    out = []
    for v in range(1 << M):
        ok = True
        for i in range(M):
            if not (v >> i) & 1 and not (v >> ((i + 1) % M)) & 1:
                ok = False
                break
        if ok:
            out.append(v)
    return np.array(sorted(out), dtype=np.int64)


class Chain:
    """Periodic golden chain of M links/anyons: basis + cached operators."""

    def __init__(self, M):
        self.M = M
        self.states = fib_states_int(M)
        self.D = self.states.size
        self.index = {int(s): i for i, s in enumerate(self.states)}
        self.cfg = np.array([[(int(s) >> i) & 1 for i in range(M)]
                             for s in self.states], dtype=np.int8)
        self._e = {}
        self._g = {}

    # --- Temperley-Lieb / braid generators (golden_full.py construction) ---

    def e(self, j):
        """TL generator at link j (sparse, real).  e^2 = phi e."""
        if j in self._e:
            return self._e[j]
        M, D = self.M, self.D
        jm, jp = (j - 1) % M, (j + 1) % M
        rows, cols, vals = [], [], []
        for col in range(D):
            s = self.cfg[col]
            a, c = s[jm], s[jp]
            if a != c:
                continue
            if a == 0:                       # vacuum neighbours: link j = tau, 1-dim
                rows.append(col); cols.append(col); vals.append(phi)
            else:                            # tau neighbours: 2-dim block phi |w><w|
                w = Fttt[:, 0]               # w_b = F[b, vac] = (1/phi, 1/sqrt(phi))
                for bnew in (0, 1):
                    t = int(self.states[col])
                    t = (t & ~(1 << j)) | (bnew << j)
                    ti = self.index.get(t)
                    if ti is not None:
                        rows.append(ti); cols.append(col)
                        vals.append(phi * w[bnew] * w[s[j]])
        E = sp.csr_matrix((vals, (rows, cols)), shape=(D, D), dtype=np.float64)
        self._e[j] = E
        return E

    def g(self, j):
        """Braid generator exchanging anyons (j, j+1): g = A + A^{-1} e_j."""
        if j not in self._g:
            E = self.e(j).astype(np.complex128)
            self._g[j] = (A_KAUFF * sp.identity(self.D, format='csr',
                                                dtype=np.complex128)
                          + (1 / A_KAUFF) * E).tocsr()
        return self._g[j]


# ------------------------------------------------- block channel projector

def _adm(a, b, c):
    """Fibonacci fusion vertex admissibility, labels in {0=vac, 1=tau}."""
    if a == 0:
        return b == c
    if b == 0:
        return a == c
    if c == 0:
        return a == b
    return True


def _fsym(a, b, c, d, e, f):
    """[F^{abc}_d]_{e f}: recoupling ((ab)_e c)_d -> (a(bc)_f)_d.
    Nontrivial 2x2 only for a=b=c=d=tau; otherwise 0/1 by admissibility."""
    if not (_adm(a, b, e) and _adm(e, c, d) and _adm(b, c, f) and _adm(a, f, d)):
        return 0.0
    if a == 1 and b == 1 and c == 1 and d == 1:
        return Fttt[e, f]
    return 1.0


def _trees(k):
    """Fusion paths c_1=tau -> ... -> c_k=vacuum with c_{j+1} in c_j x tau
    (no two consecutive vacua).  len = Fib(k-1)."""
    paths = [[1]]
    for _ in range(k - 1):
        new = []
        for p in paths:
            for c in (0, 1):
                if not (p[-1] == 0 and c == 0):
                    new.append(p + [c])
        paths = new
    return [p for p in paths if p[-1] == 0]


def block_vacuum_projector(chain, p, q):
    """P^(1) of the contiguous anyon block [p..q] (anyons p, p+1, ..., q),
    i.e. projector onto collective fusion channel = vacuum.

    F-move zip: with x0 = link (p-1) and boundary right link (q mod M),
    absorbing anyon p+j carries the pair (c_j, l_{p+j-1}) to (c_{j+1}, l_{p+j})
    with amplitude [F^{x0, c_j, tau}_{l_{p+j}}]_{l_{p+j-1}, c_{j+1}};
    the zip is unitary level by level, so  P = Z^dag Z  with Z the rows of
    the fused-basis change with root channel = vacuum.  Exact projector by
    construction.  Notes eq. (1): P_S = F_S^dag P_root F_S.
    """
    M, D = chain.M, chain.D
    k = q - p + 1
    # q > M is allowed: the block wraps through the seam (links taken mod M),
    # e.g. the defect-containing triple (N, defect, 1) on the twisted chain
    assert 1 <= p <= q and k < M, "block length must be < M"
    inner = [(p + j) % M for j in range(k - 1)]       # links p .. q-1
    lb, rb = (p - 1) % M, q % M                       # boundary links (read-only)
    inner_mask = 0
    for l in inner:
        inner_mask |= (1 << l)
    trees = _trees(k)

    groups = {}
    for i in range(D):
        key = int(chain.states[i]) & ~inner_mask
        groups.setdefault(key, []).append(i)

    rows, cols, vals = [], [], []
    for key, members in groups.items():
        x0 = (key >> lb) & 1
        y0 = (key >> rb) & 1
        # amplitude matrix A[state, tree]
        Amat = np.zeros((len(members), len(trees)))
        for a_i, i in enumerate(members):
            s = chain.cfg[i]
            ell = [s[lb]] + [s[l] for l in inner] + [s[rb]]   # l_{p-1} .. l_q
            for t_i, tr in enumerate(trees):
                amp = 1.0
                for j in range(1, k):                  # absorb anyon p+j
                    amp *= _fsym(x0, tr[j - 1], 1, ell[j + 1], ell[j], tr[j])
                    if amp == 0.0:
                        break
                Amat[a_i, t_i] = amp
        if y0 != x0:
            continue                                   # vacuum channel impossible
        P_grp = Amat @ Amat.T
        for a_i, i in enumerate(members):
            for b_i, i2 in enumerate(members):
                v = P_grp[a_i, b_i]
                if abs(v) > 1e-15:
                    rows.append(i); cols.append(i2); vals.append(v)
    return sp.csr_matrix((vals, (rows, cols)), shape=(D, D), dtype=np.float64)


# --------------------------------------------- general subset projector

def transport(chain, S, over=True, compress='right'):
    """Unitary G_S (dense) braiding the anyons of S into a contiguous block.
    Returns (G, (p, q)) with [p..q] the target block.
    over=True: moving S-anyons exchange via g (over); False: via g^dagger."""
    S = sorted(S)
    k = len(S)
    D = chain.D
    G = np.eye(D, dtype=np.complex128)

    def braid(j):
        return chain.g(j) if over else chain.g(j).getH().tocsr()

    if compress == 'right':
        targets = [S[-1] - (k - 1) + i for i in range(k)]
        for i in range(k - 2, -1, -1):                 # move right-to-left first
            for a in range(S[i], targets[i]):          # position a -> a+1 via g_a
                G = braid(a) @ G
    elif compress == 'left':
        targets = [S[0] + i for i in range(k)]
        for i in range(1, k):                          # move left-to-right first
            for a in range(S[i], targets[i], -1):      # position a -> a-1 via g_{a-1}
                G = braid(a - 1) @ G
    else:
        raise ValueError(compress)
    return np.asarray(G), (targets[0], targets[-1])


def P1_subset(chain, S, over=True, compress='right'):
    """Identity-channel projector P^(1)_S of the anyons in S (notes eq. (1)).
    Dense complex Hermitian (D, D); real iff S is contiguous.  over=True is
    the '+' chirality class; over=False gives P^- = (P^+)^*.  S within the
    fundamental domain; the arc [min S, max S] must not wrap the seam."""
    S = sorted(S)
    if all(S[i + 1] == S[i] + 1 for i in range(len(S) - 1)):   # contiguous
        return block_vacuum_projector(chain, S[0], S[-1]).toarray() \
                 .astype(np.complex128)
    G, (p, q) = transport(chain, S, over=over, compress=compress)
    Pb = block_vacuum_projector(chain, p, q)
    P = G.conj().T @ (Pb @ G)
    return 0.5 * (P + P.conj().T)


# --------------------------------------------------------- W_tau (lasso)

@njit(cache=True)
def _wtau_row(M, x_int, states_int, F):
    """Nonzero entries in column x of W_tau (verbatim from
    fib_chain/fib_chain_L22.py wtau_row, Feiguin et al. MPO)."""
    x = np.empty(M, dtype=np.int8)
    for i in range(M):
        x[i] = (x_int >> i) & 1
    free = np.empty(M, dtype=np.int64)
    n_free = 0
    for i in range(M):
        if x[i] == 1:
            free[n_free] = i
            n_free += 1
    cap = 1 << n_free
    out_idx = np.empty(cap, dtype=np.int64)
    out_val = np.empty(cap, dtype=np.float64)
    count = 0
    y = np.empty(M, dtype=np.int8)
    for mask in range(cap):
        for i in range(M):
            if x[i] == 0:
                y[i] = 1
        for j in range(n_free):
            y[free[j]] = np.int8((mask >> j) & 1)
        ok = True
        for i in range(M):
            ip1 = i + 1
            if ip1 == M:
                ip1 = 0
            if y[i] == 0 and y[ip1] == 0:
                ok = False
                break
        if not ok:
            continue
        y_int = 0
        for i in range(M):
            y_int |= np.int64(y[i]) << i
        lo, hi = 0, states_int.size
        while lo < hi:
            mid = (lo + hi) >> 1
            if states_int[mid] < y_int:
                lo = mid + 1
            else:
                hi = mid
        if lo >= states_int.size or states_int[lo] != y_int:
            continue
        val = 1.0
        for i in range(M):
            ip = i - 1
            if ip < 0:
                ip = M - 1
            if x[ip] == 1 and y[i] == 1:
                val *= F[x[i], y[ip]]
        out_idx[count] = lo
        out_val[count] = val
        count += 1
    return out_idx[:count], out_val[:count]


def build_Wtau(chain):
    """Topological symmetry / lasso W_tau on the periodic chain (dense real,
    symmetrized; asserts the fusion-algebra identity W^2 = W + 1)."""
    D = chain.D
    rows, cols, vals = [], [], []
    for c in range(D):
        idxs, vs = _wtau_row(chain.M, chain.states[c], chain.states, Fttt)
        rows.extend(idxs.tolist())
        cols.extend([c] * len(idxs))
        vals.extend(vs.tolist())
    W = sp.csr_matrix((vals, (rows, cols)), shape=(D, D)).toarray()
    asym = np.abs(W - W.T).max()
    assert asym < 1e-10, f"W_tau not symmetric: {asym:.2e}"
    W = 0.5 * (W + W.T)
    fus = np.linalg.norm(W @ W - W - np.eye(D)) / D
    assert fus < 1e-12, f"W^2 != W + 1: {fus:.2e}"
    return W


# ------------------------------------------------------------ SYK model

def all_subsets(N, k):
    """All k-subsets of physical anyons {1..N} (fundamental domain: the arc
    [min S, max S] never crosses the seam link 0 by construction)."""
    from itertools import combinations
    return [list(s) for s in combinations(range(1, N + 1), k)]


def couplings(N, k, n_subsets, seed):
    """i.i.d. Gaussian J_S, variance (k-1)!/N^(k-1) (notes eq. (2))."""
    from math import factorial
    rng = np.random.default_rng(seed)
    sd = sqrt(factorial(k - 1) / N ** (k - 1))
    return sd * rng.standard_normal(n_subsets)


def syk_hamiltonians(chain, subsets, J_list, progress=False):
    """H_r = sum_S J_S^{(r)} P^(1)_S ('+' chirality) for several coupling
    vectors at once (each P_S built once).  Dense complex Hermitian."""
    D = chain.D
    Hs = [np.zeros((D, D), dtype=np.complex128) for _ in J_list]
    for si, S in enumerate(subsets):
        P = P1_subset(chain, S)
        for r, J in enumerate(J_list):
            Hs[r] += J[si] * P
        if progress and (si + 1) % 100 == 0:
            print(f"    subset {si + 1}/{len(subsets)}")
    return [0.5 * (H + H.conj().T) for H in Hs]


# ------------------------------------------------------- tube-sector tools

def sector2_projector(W):
    """Pi_2 = (W + 1/phi) / (phi + 1/phi): kinematic, polynomial in W_tau."""
    D = W.shape[0]
    return (W + np.eye(D) / phi) / (phi + 1 / phi)


def split_sectors_13(W, H_ref):
    """Resolve sectors (1) and (3) inside the W = -1/phi block using a
    reference symmetric Hamiltonian: sector (3) = the exact energy-degenerate
    partners of sector (2); sector (1) = remainder.
    (Intrinsic construction of fib_chain/exploration/golden_Lttau.py.)
    Returns (P1, P3, max_match) with max_match the worst (2)-(3) energy match."""
    wv, wU = np.linalg.eigh(W)
    plus = wU[:, wv > 0]
    minus = wU[:, wv < 0]
    E2 = np.linalg.eigvalsh(plus.conj().T @ H_ref @ plus)
    em, vm = np.linalg.eigh(minus.conj().T @ H_ref @ minus)
    used = set()
    idx3 = []
    maxmatch = 0.0
    for E in E2:
        order = np.argsort(np.abs(em - E))
        for j in order:
            if j not in used:
                used.add(j)
                idx3.append(j)
                maxmatch = max(maxmatch, abs(em[j] - E))
                break
    idx1 = [j for j in range(len(em)) if j not in used]
    V3 = minus @ vm[:, idx3]
    V1 = minus @ vm[:, idx1]
    return V1 @ V1.conj().T, V3 @ V3.conj().T, maxmatch


def tau_cyc(chain):
    """Cyclic link shift a_j -> a_{j+1}; satisfies tau e_j tau^-1 = e_{j+1}."""
    D, M = chain.D, chain.M
    T = np.zeros((D, D))
    for col in range(D):
        s = chain.cfg[col]
        t = 0
        for b in range(M):
            t |= int(s[(b - 1) % M]) << b
        T[chain.index[t], col] = 1.0
    return T


def twisted_translation(chain, dagger=True):
    """T_tw = g_0^(dag) tau_cyc on the marked-defect chain: translates the
    physical anyons by one site, returns the defect to the seam with one
    braid (fib_chain/exploration/twisted_tube_complete.py construction)."""
    g0 = chain.g(0).toarray()
    return (g0.conj().T if dagger else g0) @ tau_cyc(chain)


def minimal_arc_rotation(S, n_phys):
    """Rotation m such that the minimal arc of the circular subset S
    (largest gap excluded) fits in the fundamental domain after shifting
    every element by +m (mod n_phys, relabelled to {1..n_phys}).
    Ties between equal largest gaps are broken INTRINSICALLY (by the
    lexicographically largest gap sequence read clockwise from the
    candidate), so the convention is rotation-covariant except for
    subsets with a nontrivial rotation symmetry, where any choice breaks
    the subset's own symmetry (documented; measure-zero fraction).
    Returns (m, S_shifted)."""
    Ss = sorted(S)
    k = len(Ss)
    gaps = [(Ss[(i + 1) % k] - Ss[i]) % n_phys for i in range(k)]
    gmax = max(gaps)
    cands = [i for i, g in enumerate(gaps) if g == gmax]
    if len(cands) > 1:
        def key(i):
            return tuple(gaps[(i + 1 + j) % k] for j in range(k - 1))
        cands.sort(key=key, reverse=True)
    imax = cands[0]
    first = Ss[(imax + 1) % k]
    m = (1 - first) % n_phys
    S_sh = sorted(((s - 1 + m) % n_phys) + 1 for s in Ss)
    return m, S_sh


def P1_subset_periodic(chain, S, n_phys, R_op=None):
    """Identity-channel projector with the PERIODIC minimal-arc convention:
    the k anyons are fused along their minimal arc, wrapped subsets handled
    by conjugation with the translation operator R_op (tau_cyc on the
    untwisted chain; twisted_translation on the marked-defect chain, which
    dresses defect crossings with the consistent chirality automatically).
    R_op moves anyon j -> j+1; P_wrapped = R^{-m} P_unwrapped(S+m) R^{m}."""
    m, S_sh = minimal_arc_rotation(S, n_phys)
    P = P1_subset(chain, S_sh)
    if m == 0:
        return P
    assert R_op is not None, "wrapped subset needs the translation operator"
    Rm = np.linalg.matrix_power(R_op, m)
    return Rm.conj().T @ P @ Rm


def _components(adj):
    n = adj.shape[0]
    labels = -np.ones(n, dtype=int)
    c = 0
    for s in range(n):
        if labels[s] >= 0:
            continue
        stack = [s]
        labels[s] = c
        while stack:
            u = stack.pop()
            for v in np.where(adj[u] | adj[:, u])[0]:
                if labels[v] < 0:
                    labels[v] = c
                    stack.append(v)
        c += 1
    return c, labels


def canonical_tube_projectors(HB0, HB1, W, dims_expect=None):
    """Invariant-block decomposition of the algebra generated by two generic
    symmetric Hamiltonians: X commuting with the (nondegenerate) HB0 is
    diagonal in its eigenbasis; [X, HB1] = 0 then forces X constant on the
    connected components of the graph with edges where <a|HB1|b> != 0.
    Returns ([{dim, P, w}], n_components); the components are the minimal
    invariant blocks (= sector projectors when the algebra is the full
    symmetric algebra)."""
    ev0, V0 = np.linalg.eigh(HB0)
    H1t = V0.conj().T @ HB1 @ V0
    adj = np.abs(H1t) > 1e-8 * np.abs(H1t).max()
    n_comp, labels = _components(adj)
    Pis = []
    for c in range(n_comp):
        idx = np.where(labels == c)[0]
        Vc = V0[:, idx]
        P = Vc @ Vc.conj().T
        w_val = np.trace(W @ P).real / len(idx)
        Pis.append({'dim': len(idx), 'P': P, 'w': w_val})
    return Pis, n_comp


def r_stat(E, drop=0.05):
    """Mean adjacent-gap ratio <r> (GOE .5307, GUE .5996, GSE .6744, Poi .3863).
    Degenerate levels must be removed by the caller (sector-resolved input)."""
    E = np.sort(E)
    n = len(E)
    lo, hi = int(drop * n), int((1 - drop) * n)
    E = E[lo:hi]
    s = np.diff(E)
    r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return r


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a
