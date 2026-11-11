# Data

## `data/fib/fib_L26_lam1.25_spectra.npz`

Energy eigenvalues of the perturbed Fibonacci golden chain of Sec. 2.5.2 of the paper. The chain has $L = 26$ anyons and coupling $\lambda = 1.25$, and the eigenvalues are resolved into the symmetry blocks plotted in Fig. 7 (`numerics_fib`). The file is 3.7 MB and holds 102 one-dimensional float64 arrays, each sorted in increasing order. It is required for the default workflow and needs no download.

| Key | Block | Blocks stored | Levels per block |
|---|---|---|---|
| `U\|k<k>W+\|E` | untwisted, $\mathcal L^{\tau,\tau}_{\mathbf 1\leftarrow\mathbf 1}$ eigenvalue $\varphi$, center label $(\mathbf 1,\mathbf 1)$ ($W_+$) | $k = 1,\dots,12$ | 2871 to 2900 |
| `U\|k<k>W-\|E` | untwisted, eigenvalue $-1/\varphi$, center label $(\tau,\tau)$ ($W_-$) | $k = 1,\dots,12$ | 7549 to 7560 |
| `T\|s1m<m>\|E` | $\tau$-twisted, $(\tau,\tau)$ | 26 values of $m$ | 7549 to 7561 |
| `T\|s2m<m>\|E` | $\tau$-twisted, $(\tau,\mathbf 1)$, spin $e^{4\pi i/5}$ | 26 values of $m$ | 4660 to 4678 |
| `T\|s3m<m>\|E` | $\tau$-twisted, $(\mathbf 1,\tau)$, spin $e^{-4\pi i/5}$ | 26 values of $m$ | 4660 to 4678 |

Here $k$ is the lattice momentum index of the untwisted chain, $2\pi k/L$. The untwisted momenta stored are one from each $\pm k$ pair. The self-conjugate momenta $k = 0$ and $k = 13$, where reflection acts as an extra unitary symmetry, are not stored (Sec. 2.5.2).

In the twisted chain, $m$ labels the eigenvalue $e^{2\pi i m/(5L)}$ of the twisted translation, which generates $\mathbb Z_{5L}$, with $m$ taken mod $5L = 130$. Its $L$-th power is the $2\pi$ rotation and gives the topological spin $e^{2\pi i m/5}$, so $m \equiv 0$, $2$ and $3 \pmod 5$ for $(\tau,\tau)$, $(\tau,\mathbf 1)$ and $(\mathbf 1,\tau)$. The untwisted block $W_-$ at momentum $k$ has the same spectrum as the twisted $(\tau,\tau)$ block at $m = 5k$. The $(\tau,\mathbf 1)$ block at $m$ has the same spectrum as the $(\mathbf 1,\tau)$ block at $130 - m$. `code/figure-reproduction/check_fib_isospectrality.py` uses these matchings.

The full twisted Hilbert space is stored, so the twisted totals are exact. Over the 26 momenta, $(\tau,\tau)$ adds up to $F_{27} = 196418$ levels, and $(\tau,\mathbf 1)$ and $(\mathbf 1,\tau)$ each add up to $F_{26} = 121393$.

**Provenance.** The eigenvalues were computed in July 2026 on the Fir cluster of the Digital Research Alliance of Canada. The programs were `code/fib/fib_eth_lowmem.py` (untwisted chain) and `code/fib/fib_eth_twisted_lowmem.py` (twisted chain, one job per block). The original output files also contain eigenvectors and operator matrix elements. Those are used in the companion paper on categorical eigenstate thermalization and are not included here. `code/fib/extract_fib_L26_spectra.py` copied every eigenvalue array from those files into this one, with no other change. The edge trimming and unfolding of Fig. 7 happen later, in the figure script.

**Used by:**
- `code/figure-reproduction/fig07_fib_sectors_L26.py` (Fig. 7 and the five values of $\langle r\rangle$ in Sec. 2.5.2)
- `code/figure-reproduction/check_fib_isospectrality.py` (the level-by-level coincidences stated in Sec. 2.5.2).

**Regenerating it.** A fresh computation at $L = 26$ needs cluster-scale memory and is not part of the default workflow. `code/figure-reproduction/check_fib_generators_smallL.py` runs the same two programs at $L = 14$ in seconds. `code/figure-reproduction/README.md` gives the details.

## `data/s3/`

Energy eigenvalues of the disordered $S_3$ qutrit chain of Sec. 2.4.3, with the Hamiltonian `eq:S3ham` and the disorder `eq:S3disorder`, at $L = 8$. Each file holds 120 disorder realizations, and in each realization the spectra of all eight $D(S_3)$ sectors, with the same couplings for the three twists $e$, $\sigma$ and $r$. Both files are required for the default workflow and need no download.

| File | Ensemble | Seeds, realization $r = 0,\dots,119$ | Used for | Size |
|---|---|---|---|---|
| `S3_rep_category_L8_spectra.npz` | June 2026 | numpy `default_rng(20248601 + r)` | Table 1, Figs. 3 and 4 | 17 MB |
| `ensemble_L8_n120.npz` | November 2026 | numpy `default_rng(20261103 + r)` | Figs. 5 and 6 | 16 MB |

The sectors $([e],1)$, $([e],\mathrm{sgn})$, $([e],2)$, $([t],+)$, $([t],-)$, $([c],1)$, $([c],\omega)$ and $([c],\omega^2)$ have 1094, 1093, 2187, 3281, 3280, 2187, 2187 and 2187 levels.

**Formats.** The June file holds 960 one-dimensional float64 arrays, one per sector and realization, named `<code>__<r>` with codes `e_1`, `e_sgn`, `e_2`, `t_+`, `t_-`, `c_1`, `c_w` and `c_w2`. It holds spectra only. The November file holds one array `spec_<sector>` of shape (120, levels) per sector, with sector names `e1`, `esgn`, `e2`, `t+`, `t-`, `c1`, `cw` and `cw2`. It also holds the couplings of every realization (`coup_J`, `coup_f`, `coup_J2`, `coup_J3`), the seed base, and the numpy version and platform that made it. `load_npz` in `code/s3/python/run_paper_numerics.py` reads both formats.

**Two ensembles.** As Sec. 2.4.3 states, the paper uses two independent sets of 120 realizations: Table 1 and Figs. 3 and 4 come from the June ensemble, and Figs. 5 and 6 from the November one. (arXiv v1 described a single set.) They are independent draws from the same disorder distribution and give the same picture: the gap ratios of the November ensemble agree with Table 1 within 1.5 combined standard errors, and in both ensembles the off-diagonal $C_{AB}$ stay below 0.03. `code/figure-reproduction/s3_from_stored_ensembles.py` prints these comparisons.

**Provenance.** Both files were computed by Sridip Pal. The November file was made with `code/s3/python/run_paper_numerics.py` on an Apple M4 with Python 3.9.6 and numpy 2.0.2, as recorded in the file. The June file was made in June 2026 by an earlier script that is not included, and it records no seeds. That run's seed rule, 20240601 + 1000 $L$ + $r$, is seed base 20248601 at $L = 8$, and with it the included code regenerates the June spectra (see below).

**Used by:**
- `code/figure-reproduction/s3_from_stored_ensembles.py` (Table 1, Figs. 3, 5 and 6, and the quick checks of Sec. 2.4.3)
- `code/figure-reproduction/fig04_s3_sff_goe_binned.py` (Fig. 4)
- `code/figure-reproduction/check_s3_regeneration.py` (comparison with freshly regenerated spectra).

**Regenerating them.** `code/figure-reproduction/check_s3_regeneration.py june` (or `nov`) recomputes an ensemble from its seeds in about 9 minutes on the tested Mac and compares it level by level with the stored file. Regeneration relies on numpy's random stream, which numpy does not formally guarantee across versions, so the stored files are the permanent record.
