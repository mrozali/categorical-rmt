"""
Spectral diagnostics of "Categorical RMT", Section 2.4.3, paragraph "Diagnostics":
  * unfolding: degree-8 polynomial fit to the level staircase N(E) = #{n : E_n <= E}; the unfolded
    levels eps_n = poly(E_n) have unit mean density;
  * every diagnostic uses the central 70% of each spectrum (by level index);
  * P(s): distribution of unfolded gaps s_n = eps_{n+1} - eps_n, against the Wigner surmises;
  * adjacent-gap ratio r_n = min(d_n, d_{n+1}) / max(d_n, d_{n+1}) from the RAW gaps d_n = E_{n+1} - E_n
    (unfolding independent), against <r> = 0.38629 (Poisson), 0.5307 (GOE), 0.5996 (GUE);
  * connected SFF, eq. (2.25):  K_c(t)/N = (<|Z|^2> - |<Z>|^2) / <N>,  Z(t) = sum_n exp(-2 pi i eps_n t),
    with time in units of the Heisenberg time (t_H = 1 for unit mean density);
  * cross SFF, eq. (2.26): Re(<Z_A Z_B^*> - <Z_A><Z_B>^*) / sqrt(<N_A><N_B>);
  * correlation matrix, eq. (2.27): with dZ_X the connected part of N_X^{-1/2} Z_X,
    C_AB = <Re dZ_A dZ_B^*> / sqrt(<|dZ_A|^2><|dZ_B|^2>), averaged over realisations and a band of times.
RMT references for the SFF (unit density, t in units of t_H; standard results, see e.g. M. L. Mehta,
Random Matrices, 3rd ed., Elsevier 2004):
    GUE: K/N = t (t<1), 1 (t>=1);   GOE: K/N = 2t - t ln(1+2t) (t<1),  2 - t ln((2t+1)/(2t-1)) (t>1).
"""
from __future__ import annotations
import numpy as np

R_POISSON, R_GOE, R_GUE = 2 * np.log(2) - 1, 0.5307, 0.5996
R_REF = {"Poisson": R_POISSON, "GOE": R_GOE, "GUE": R_GUE}

def central_window(n, frac=0.7):
    lo = int(round(n * (1 - frac) / 2)); return lo, n - lo

def unfold(E, deg=8, frac=0.7):
    """unfolded levels of the central window, from a degree-`deg` polynomial fit of the full staircase."""
    E = np.sort(E); n = len(E)
    N = np.arange(1, n + 1)
    # fit on a rescaled variable for conditioning
    x = (E - E.mean()) / E.std()
    c = np.polyfit(x, N, deg)
    lo, hi = central_window(n, frac)
    return np.polyval(c, x[lo:hi])

def gap_ratios(E, frac=0.7):
    """r_n from raw gaps in the central window."""
    E = np.sort(E); lo, hi = central_window(len(E), frac)
    d = np.diff(E[lo:hi])
    return np.minimum(d[:-1], d[1:]) / np.maximum(d[:-1], d[1:])

def rbar_ensemble(spectra, frac=0.7):
    """mean gap ratio over an ensemble (nreal, dim): returns (mean, sem over realisations, sem pooled, n_r)."""
    per_real = []; pooled = []
    for E in spectra:
        rr = gap_ratios(E, frac); per_real.append(rr.mean()); pooled.append(rr)
    per_real = np.array(per_real); pooled = np.concatenate(pooled)
    return per_real.mean(), per_real.std(ddof=1) / np.sqrt(len(per_real)), pooled.std(ddof=1) / np.sqrt(len(pooled)), len(pooled)

def surmise(s, cls):
    if cls == "Poisson": return np.exp(-s)
    if cls == "GOE": return (np.pi / 2) * s * np.exp(-np.pi * s ** 2 / 4)
    if cls == "GUE": return (32 / np.pi ** 2) * s ** 2 * np.exp(-4 * s ** 2 / np.pi)
    raise ValueError(cls)

def sff_rmt(t, cls):
    t = np.asarray(t, dtype=float)
    if cls == "GUE":
        return np.where(t < 1, t, 1.0)
    if cls == "GOE":
        small = 2 * t - t * np.log1p(2 * t)
        with np.errstate(divide="ignore", invalid="ignore"):
            large = 2 - t * np.log((2 * t + 1) / (2 * t - 1))
        return np.where(t < 1, small, large)
    raise ValueError(cls)

def amplitudes(spectra, ts, frac=0.7, deg=8):
    """Z_r(t) = sum_n exp(-2 pi i eps_n t) for every realisation r (rows) and time (columns), and N_r."""
    Z = np.empty((len(spectra), len(ts)), dtype=complex); Ns = np.empty(len(spectra))
    for r, E in enumerate(spectra):
        eps = unfold(E, deg, frac); Ns[r] = len(eps)
        Z[r] = np.exp(-2j * np.pi * np.outer(ts, eps)).sum(axis=1)
    return Z, Ns

def sff_connected(Z, Ns):
    """K_c(t)/N of eq. (2.25)."""
    return (np.mean(np.abs(Z) ** 2, axis=0) - np.abs(Z.mean(axis=0)) ** 2) / Ns.mean()

def cross_sff(ZA, NA, ZB, NB):
    """K_c^{AB}(t)/N of eq. (2.26)."""
    return np.real(np.mean(ZA * np.conj(ZB), axis=0) - ZA.mean(axis=0) * np.conj(ZB.mean(axis=0))) / np.sqrt(NA.mean() * NB.mean())

def corr_matrix(Zs, Ns, ts, band):
    """C_AB of eq. (2.27) on the time band [band[0], band[1]] (inclusive), averaged over
    realisations and times.  Zs: dict sector -> (nreal, nt) amplitudes; Ns: dict sector -> N per realisation."""
    keep = (ts >= band[0]) & (ts <= band[1])
    dZ = {}
    for k in Zs:
        z = Zs[k][:, keep] / np.sqrt(Ns[k])[:, None]
        dZ[k] = z - z.mean(axis=0, keepdims=True)
    keys = list(Zs); C = np.zeros((len(keys), len(keys)))
    for i, a in enumerate(keys):
        for j, b in enumerate(keys):
            num = np.mean(np.real(dZ[a] * np.conj(dZ[b])))
            C[i, j] = num / np.sqrt(np.mean(np.abs(dZ[a]) ** 2) * np.mean(np.abs(dZ[b]) ** 2))
    return keys, C, int(keep.sum())
