"""
S3-symmetric qutrit chain with twisted boundary conditions of "Categorical RMT" (Cui, Jensen, Pal, Rozali),
Section 2.4: the twisted Hamiltonian of eq. (2.23), the disorder of eq. (2.24), and the bases of the eight
tube algebra sectors ([h],R) of D(S3) defined by the projectors of eqs. (2.17)-(2.18).

Conventions (all of them taken from the manuscript):
  * sites i = 0, ..., L-1 carry a qutrit |s_i>, s_i in {0,1,2}; the many-body basis state with labels
    (s_0, ..., s_{L-1}) is the integer n = sum_i s_i 3^i  (site i is digit i of n in base 3);
  * S3 acts by permuting the labels, U(g)|s_0,...,s_{L-1}> = |g(s_0),...,g(s_{L-1})>, with generators
    sigma = (1 2) and r = (0 -> 1 -> 2 -> 0);
  * the twist h enters only through the wrap-around bond J_{L-1} delta_{s_{L-1}, h(s_0)}; the J2 and J3
    sums stop at i = L-3 (no wrap);
  * the transverse field is f_i sum_{a != b} |a><b|_i, i.e. +f_i on every off-diagonal on-site element;
  * sectors: for h = e the three S3 irreps, for h = sigma the two signs of U_sigma(sigma), for h = r the
    three eigenvalues 1, omega, omega^2 of U_r(r) with omega = exp(2 pi i/3).

Every sector basis is built explicitly from orbits of the centraliser on basis states (an independent
route to the projector ranges), and the two-dimensional irrep is handled by two different routes,
which must agree.  Run  python3 s3_chain.py  to execute the self-tests (well under a second at L = 5).
"""
from __future__ import annotations
import itertools, sys, time
import numpy as np
import scipy.sparse as sp

# ----------------------------------------------------------------------------------------------
# the group S3 as permutations of the labels {0,1,2}; an element is the tuple (g(0), g(1), g(2))
# ----------------------------------------------------------------------------------------------
E = (0, 1, 2)
SIGMA = (0, 2, 1)          # transposition (1 2)
R = (1, 2, 0)              # three-cycle 0 -> 1 -> 2 -> 0

def mul(g, h):
    """(g h)(s) = g(h(s)): first h, then g."""
    return tuple(g[h[s]] for s in range(3))

def inv(g):
    out = [0, 0, 0]
    for s in range(3):
        out[g[s]] = s
    return tuple(out)

GROUP = sorted({g for g in itertools.permutations(range(3))})
NAME = {E: "e", SIGMA: "s", mul(R, SIGMA): "rs", mul(mul(R, R), SIGMA): "r2s", R: "r", mul(R, R): "r2"}
assert len(GROUP) == 6 and set(NAME) == set(GROUP)

def sign(g):
    """sign of the permutation g (character of the sign representation)."""
    return 1 if sum(1 for a in range(3) for b in range(a + 1, 3) if g[a] > g[b]) % 2 == 0 else -1

def chi2(g):
    """character of the two-dimensional irrep = (#fixed labels) - 1."""
    return sum(1 for s in range(3) if g[s] == s) - 1

def centraliser(h):
    return [g for g in GROUP if mul(g, h) == mul(h, g)]

# real orthogonal matrices of the two-dimensional irrep: the permutation representation of S3 on R^3
# restricted to the plane orthogonal to (1,1,1), in the orthonormal basis e1 = (1,-1,0)/sqrt2,
# e2 = (1,1,-2)/sqrt6.  (P(g) v)_s = v_{g^{-1}(s)}.
_E1 = np.array([1, -1, 0]) / np.sqrt(2.0)
_E2 = np.array([1, 1, -2]) / np.sqrt(6.0)
def D2(g):
    gi = inv(g)
    P = np.zeros((3, 3))
    for s in range(3):
        P[s, gi[s]] = 1.0          # (P v)_s = v_{g^{-1}(s)}
    B = np.stack([_E1, _E2], axis=1)
    return B.T @ P @ B

# ----------------------------------------------------------------------------------------------
# many-body basis and the group action
# ----------------------------------------------------------------------------------------------
def digits(L):
    """array of shape (3^L, L): digits[n, i] = s_i of basis state n."""
    n = np.arange(3 ** L)
    return np.stack([(n // 3 ** i) % 3 for i in range(L)], axis=1).astype(np.int8)

def perm_index(L, g, dig=None):
    """index of U(g)|n> for every n: apply g to every label."""
    dig = digits(L) if dig is None else dig
    gd = np.asarray(g, dtype=np.int64)[dig]
    return (gd * (3 ** np.arange(L))[None, :]).sum(axis=1)

def U(L, g, dig=None):
    """permutation matrix of U(g) as a sparse matrix: U|n> = |g.n>, i.e. U[g.n, n] = 1."""
    N = 3 ** L
    tgt = perm_index(L, g, dig)
    return sp.csr_matrix((np.ones(N), (tgt, np.arange(N))), shape=(N, N))

# ----------------------------------------------------------------------------------------------
# the twisted Hamiltonian, eq. (2.23), and the disorder of eq. (2.24)
# ----------------------------------------------------------------------------------------------
COUPLINGS = {"J": (1.00, 0.15), "f": (0.82, 0.12), "J2": (0.32, 0.08), "J3": (0.20, 0.06)}

def draw_couplings(L, rng):
    """one disorder realisation of eq. (2.24): mean + width * standard normal, per site and term."""
    m = COUPLINGS
    J = m["J"][0] + m["J"][1] * rng.standard_normal(L)
    f = m["f"][0] + m["f"][1] * rng.standard_normal(L)
    J2 = m["J2"][0] + m["J2"][1] * rng.standard_normal(L - 2)
    J3 = m["J3"][0] + m["J3"][1] * rng.standard_normal(L - 2)
    return J, f, J2, J3

def hamiltonian(L, h, J, f, J2, J3, dig=None):
    """H_h of eq. (2.23) as a real symmetric sparse matrix on the 3^L-dimensional space."""
    dig = digits(L) if dig is None else dig
    N = 3 ** L
    s = dig.astype(np.int64)
    hmap = np.asarray(h, dtype=np.int64)
    # diagonal part: Potts bonds (with the twisted wrap bond), next-nearest-neighbour, three-site
    diag = np.zeros(N)
    for i in range(L - 1):
        diag += J[i] * (s[:, i] == s[:, i + 1])
    diag += J[L - 1] * (s[:, L - 1] == hmap[s[:, 0]])          # delta_{s_{L-1}, h(s_0)}
    for i in range(L - 2):
        diag += J2[i] * (s[:, i] == s[:, i + 2])
        diag += J3[i] * ((s[:, i] == s[:, i + 1]) & (s[:, i + 1] == s[:, i + 2]))
    rows = [np.arange(N)]; cols = [np.arange(N)]; vals = [diag]
    # transverse field: f_i sum_{a != b} |a><b|_i  (every change of one site label, amplitude f_i)
    n = np.arange(N)
    for i in range(L):
        for b in range(3):
            mask = s[:, i] != b
            m = n[mask]
            target = m + (b - s[mask, i]) * 3 ** i
            rows.append(target); cols.append(m); vals.append(np.full(m.size, f[i]))
    H = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return H

# ----------------------------------------------------------------------------------------------
# sector bases, built from orbits of the centraliser on basis states
# ----------------------------------------------------------------------------------------------
SECTORS = ["e1", "esgn", "e2", "t+", "t-", "c1", "cw", "cw2"]
SECTOR_TWIST = {"e1": E, "esgn": E, "e2": E, "t+": SIGMA, "t-": SIGMA, "c1": R, "cw": R, "cw2": R}
SECTOR_CLASS = {"e1": "GOE", "esgn": "GOE", "e2": "GOE", "t+": "GOE", "t-": "GOE", "c1": "GOE",
                "cw": "GUE", "cw2": "GUE"}
OMEGA = np.exp(2j * np.pi / 3)

def _orbits(L, elements, dig):
    """orbits of the subgroup `elements` (the complete list of its elements) on basis states.  Returns a list
    of pairs (n0, members): n0 is the smallest state of the orbit and members = [g.n0 for g in elements], in
    the order of `elements` and WITH repetitions whenever n0 has a nontrivial stabilizer.  The repetitions
    matter: summing chi(g) over members gives the unnormalized projector image P^R|n0> of eqs. (2.17)-(2.18)."""
    N = 3 ** L
    maps = {g: perm_index(L, g, dig) for g in elements}
    seen = np.zeros(N, dtype=bool); orbits = []
    for n0 in range(N):
        if seen[n0]:
            continue
        members = [int(maps[g][n0]) for g in elements]
        uniq = list(dict.fromkeys(members))
        for m in uniq:
            seen[m] = True
        orbits.append((n0, members))
    return orbits

def sector_basis(L, sector, dig=None, route="row"):
    """orthonormal basis (columns of a sparse 3^L x d matrix) of the multiplicity space of `sector` inside
    H_h, h = SECTOR_TWIST[sector].  Every column is the normalized projector image P|n0> of one orbit:
      ([e],1), ([e],sgn):  P^1_e|n0>, P^sgn_e|n0> of eq. (2.17), i.e. sum_g chi(g) U(g)|n0>;
      ([e],2):             P^2_e of eq. (2.17) projects onto C^2 (x) M_([e],2), eq. (2.16), so one copy of C^2
                           gives a block isospectral to H_([e],2).  route="row" (default) uses the matrix-element
                           projector (2/6) sum_g D_11(g) U(g), which selects the first basis vector of C^2;
                           route="branch" uses P^2_e (1 + U(sigma))/2, the U(sigma) = +1 vector, as the
                           Mathematica code does.  The self-test checks that the two routes agree;
      ([t],+), ([t],-):    P^+-_sigma|n0> = (|n0> +- U(sigma)|n0>)/2 of eq. (2.18), normalized;
      ([c],chi):           sqrt3 P^chi_r|n0> of eq. (2.18) = (1/sqrt3) sum_k chi(r)^(-k) U(r)^k |n0>."""
    dig = digits(L) if dig is None else dig
    N = 3 ** L
    cols_idx, cols_val, col = [], [], 0
    def add(vec_idx, vec_val):
        nonlocal col
        cols_idx.append((vec_idx, np.full(len(vec_idx), col))); cols_val.append(vec_val); col += 1
    h = SECTOR_TWIST[sector]
    if h == E:
        orbs = _orbits(L, GROUP, dig)
        if sector in ("e1", "esgn"):
            char = (lambda g: 1) if sector == "e1" else sign
            for n0, members in orbs:
                v = {}
                for g, m in zip(GROUP, members):
                    v[m] = v.get(m, 0.0) + char(g)
                idx = np.array(list(v)); val = np.array([v[k] for k in idx], dtype=float)
                nrm = np.linalg.norm(val)
                # the sign vector vanishes on the size-3 orbit {|aa...a>}, whose states are fixed by a
                # transposition, so that orbit gives no ([e],sgn) vector; the cut detects this exact zero
                if nrm > 1e-9:
                    add(idx, val / nrm)
        else:  # the two-dimensional irrep
            for n0, members in orbs:
                uniq = list(dict.fromkeys(members)); pos = {m: k for k, m in enumerate(uniq)}
                # matrices of U(g) restricted to the orbit, in the basis `uniq`
                Ug = {}
                for g in GROUP:
                    M = np.zeros((len(uniq), len(uniq)))
                    for m in uniq:
                        M[pos[perm_from(L, g, m, dig)], pos[m]] = 1.0
                    Ug[g] = M
                if route == "row":
                    P = (2.0 / 6.0) * sum(D2(g)[0, 0] * Ug[g] for g in GROUP)
                else:
                    P = (2.0 / 6.0) * sum(chi2(g) * Ug[g] for g in GROUP)
                    P = P @ (0.5 * (np.eye(len(uniq)) + Ug[SIGMA]))
                # orthonormal basis of the range of P (an orthogonal projector in both routes, since P^2_e is
                # central and commutes with U(sigma))
                w, V = np.linalg.eigh(P @ P.T)
                keep = V[:, w > 1e-9]
                Q, _ = np.linalg.qr(keep)
                for k in range(Q.shape[1]):
                    add(np.array(uniq), Q[:, k])
    elif h == SIGMA:
        orbs = _orbits(L, [E, SIGMA], dig); sgn = +1 if sector == "t+" else -1
        for n0, (m0, m1) in orbs:
            if m0 == m1:
                if sgn == +1:
                    add(np.array([m0]), np.array([1.0]))
            else:
                add(np.array([m0, m1]), np.array([1.0, sgn]) / np.sqrt(2))
    else:  # h == R, eigenvalues of U_r(r): 1, omega, omega^2
        lam = {"c1": 1.0, "cw": OMEGA, "cw2": OMEGA ** 2}[sector]
        orbs = _orbits(L, [E, R, mul(R, R)], dig)
        for n0, members in orbs:
            assert len(set(members)) == 3
            # v = (1/sqrt3) sum_k lam^{-k} U(r)^k |n0>  satisfies U(r) v = lam v
            add(np.array(members), np.array([lam ** (-k) for k in range(3)]) / np.sqrt(3))
    idx = np.concatenate([c[0] for c in cols_idx]); jdx = np.concatenate([c[1] for c in cols_idx])
    val = np.concatenate(cols_val)
    dtype = complex if np.iscomplexobj(val) else float
    return sp.csc_matrix((val.astype(dtype), (idx, jdx)), shape=(N, col))

def perm_from(L, g, n, dig):
    """index of U(g)|n> for a single basis state n."""
    s = dig[n].astype(np.int64)
    return int((np.asarray(g)[s] * 3 ** np.arange(L)).sum())

def block(H, B):
    """B^dagger H B as a dense matrix."""
    M = B.conj().T @ (H @ B)
    return np.asarray(M.todense()) if sp.issparse(M) else np.asarray(M)

def sector_dimension_formula(L, sector):
    """m_rho = (1/|C(h)|) sum_{g in C(h)} chi_rho(g)^* Fix(g)^L with Fix(g) = #fixed labels."""
    h = SECTOR_TWIST[sector]; C = centraliser(h)
    fix = lambda g: sum(1 for s in range(3) if g[s] == s)
    if sector == "e1": ch = lambda g: 1
    elif sector == "esgn": ch = sign
    elif sector == "e2": ch = chi2
    elif sector == "t+": ch = lambda g: 1
    elif sector == "t-": ch = lambda g: 1 if g == E else -1
    else:
        lam = {"c1": 1.0, "cw": OMEGA, "cw2": OMEGA ** 2}[sector]
        ch = lambda g: {E: 1.0, R: lam, mul(R, R): lam ** 2}[g]
    val = sum(np.conj(ch(g)) * fix(g) ** L for g in C) / len(C)
    return int(round(val.real))

# ----------------------------------------------------------------------------------------------
# self-tests
# ----------------------------------------------------------------------------------------------
def selftest(L=5, seed=1):
    rng = np.random.default_rng(seed)
    dig = digits(L); N = 3 ** L
    J, f, J2, J3 = draw_couplings(L, rng)
    t0 = time.time()
    print(f"self-test at L={L} (N={N})")
    # 1. group structure and representation
    for g in GROUP:
        for g2 in GROUP:
            assert (U(L, g, dig) @ U(L, g2, dig) != U(L, mul(g, g2), dig)).nnz == 0
    assert all(np.allclose(D2(g) @ D2(g2), D2(mul(g, g2))) for g in GROUP for g2 in GROUP)
    assert all(np.allclose(D2(g) @ D2(g).T, np.eye(2)) for g in GROUP)
    assert all(abs(np.trace(D2(g)) - chi2(g)) < 1e-12 for g in GROUP)
    print("  [ok] U(g)U(g') = U(gg'); D2 is a real orthogonal representation with character chi2")
    # 2. intertwiners and centraliser commutation, U(g) H_h U(g)^-1 = H_{g h g^-1}
    Hs = {h: hamiltonian(L, h, J, f, J2, J3, dig) for h in GROUP}
    for h in GROUP:
        assert (Hs[h] - Hs[h].T).count_nonzero() == 0
        for g in GROUP:
            lhs = U(L, g, dig) @ Hs[h] @ U(L, inv(g), dig)
            assert abs(lhs - Hs[mul(mul(g, h), inv(g))]).max() < 1e-12
            comm = abs(Hs[h] @ U(L, g, dig) - U(L, g, dig) @ Hs[h]).max()
            assert (comm < 1e-12) == (g in centraliser(h)), (NAME[h], NAME[g], comm)
    print("  [ok] H_h real symmetric; U(g) H_h U(g)^-1 = H_{ghg^-1}; [H_h, U(g)] = 0 iff g in C(h)")
    # 3. traces: Tr H_h = (N/3)(sum J + sum J2) + (N/9) sum J3 for every twist
    for h in GROUP:
        assert abs(Hs[h].diagonal().sum() - (N / 3 * (J.sum() + J2.sum()) + N / 9 * J3.sum())) < 1e-9
    print("  [ok] Tr H_h = (N/3)(sum J_i + sum J2_i) + (N/9) sum J3_i, independent of the twist")
    # 4. sector bases: orthonormal, dimensions, invariance, recombination of spectra
    spec = {}
    for sec in SECTORS:
        B = sector_basis(L, sec, dig)
        G = block(sp.identity(N), B)
        assert np.allclose(G, np.eye(B.shape[1]), atol=1e-12), sec
        assert B.shape[1] == sector_dimension_formula(L, sec), (sec, B.shape[1], sector_dimension_formula(L, sec))
        h = SECTOR_TWIST[sec]
        Hb = block(Hs[h], B)
        # the basis spans an invariant subspace: H B = B (B^dag H B)
        assert abs(Hs[h] @ B - B @ Hb).max() < 1e-10, sec
        spec[sec] = np.linalg.eigvalsh(Hb)
    print("  [ok] sector bases orthonormal and invariant; dimensions:",
          {s: int(sector_dimension_formula(L, s)) for s in SECTORS})
    full = {h: np.linalg.eigvalsh(Hs[h].toarray()) for h in (E, SIGMA, R)}
    rec_e = np.sort(np.concatenate([spec["e1"], spec["esgn"], spec["e2"], spec["e2"]]))
    rec_t = np.sort(np.concatenate([spec["t+"], spec["t-"]]))
    rec_c = np.sort(np.concatenate([spec["c1"], spec["cw"], spec["cw2"]]))
    assert np.abs(rec_e - full[E]).max() < 1e-9 and np.abs(rec_t - full[SIGMA]).max() < 1e-9
    assert np.abs(rec_c - full[R]).max() < 1e-9
    print("  [ok] sector spectra recombine into Spec H_e (with ([e],2) twice), Spec H_sigma, Spec H_r")
    # conjugate twists are isospectral
    for h in GROUP:
        for g in GROUP:
            k = mul(mul(g, h), inv(g))
            assert np.abs(np.linalg.eigvalsh(Hs[h].toarray()) - np.linalg.eigvalsh(Hs[k].toarray())).max() < 1e-9
    print("  [ok] Spec H_h = Spec H_{ghg^-1} for all conjugate twists")
    # 5. the two routes to the ([e],2) block agree, and ([c],omega^2) is the complex conjugate of ([c],omega)
    B_row = sector_basis(L, "e2", dig, "row"); B_br = sector_basis(L, "e2", dig, "branch")
    assert np.abs(np.linalg.eigvalsh(block(Hs[E], B_row)) - np.linalg.eigvalsh(block(Hs[E], B_br))).max() < 1e-9
    Bw, Bw2 = sector_basis(L, "cw", dig), sector_basis(L, "cw2", dig)
    Hw, Hw2 = block(Hs[R], Bw), block(Hs[R], Bw2)
    assert np.allclose(Hw2, Hw.conj(), atol=1e-12) and np.abs(spec["cw"] - spec["cw2"]).max() < 1e-10
    assert np.abs(Hw.imag).max() > 1e-3, "the ([c],omega) block should be genuinely complex"
    print("  [ok] ([e],2): row-projector and Z2-branching routes isospectral; H_{cw2} = conj(H_cw), equal spectra")
    # 6. time reversal Theta, eqs. (2.20)-(2.21): complex conjugation maps the U_r(r) = omega basis onto the
    #    omega^2 basis, so H_([c],omega^2) = conj(H_([c],omega)) and the two spectra coincide, eq. (2.22)
    assert abs(Bw.conj() - Bw2).max() < 1e-12
    print("  [ok] Theta (complex conjugation) maps the ([c],omega) basis to the ([c],omega^2) basis")
    print(f"self-test passed in {time.time()-t0:.1f}s")

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    selftest(L)
