# S₃ qutrit chain · Python code

**Numerics for Section 2.4.3 of "Categorical RMT"** (Cui, Jensen, Pal, Rozali)

This folder reproduces every number and figure that the manuscript reports for the S₃-symmetric qutrit chain.
The model is the twisted Hamiltonian (2.23) with the disorder (2.24), resolved into the eight tube algebra
sectors ([h],R) of D(S₃). An independent Mathematica twin lives in `../S3_mathematica_paper/`, and the two codes
are cross-checked against each other.

**Contents** · [Quick start](#quick-start) · [Files](#files) · [The three ensembles](#the-three-ensembles) ·
[What reproduces what](#what-reproduces-what) · [Results](#results) · [Validation](#validation) ·
[Usage in detail](#usage-in-detail) · [Notes on the figures](#notes-on-the-figures) ·
[Packaging](#packaging-for-a-supplement)

---

## Quick start

Requires Python ≥ 3.9 with numpy, scipy and matplotlib. Every command works from any directory.

```bash
python3 s3_chain.py                                                       # self-test               0.2 s
python3 run_paper_numerics.py --seed-base 20248601                        # Table 1, Figs. 3-4      13 min
python3 run_paper_numerics.py                                             # Figs. 5-6               18 min
python3 run_paper_numerics.py --from-npz S3_rep_category_L8_spectra.npz   # Table 1 from June file  4 s
```

Each run writes Table 1, Figures 3–6, a summary of the results and a log into its own `results*/` folder. Every option is explained in [`DOCUMENTATION.md`](DOCUMENTATION.md), also available as [`DOCUMENTATION.pdf`](DOCUMENTATION.pdf).

---

## Files

| File | Contents |
|---|---|
| `DOCUMENTATION.md` | the user guide: every option, worked examples, reference |
| `DOCUMENTATION.tex`, `.pdf` | the same guide as a LaTeX document and its PDF |
| `s3_chain.py` | the model: S₃ action, Hamiltonian (2.23), disorder (2.24), sector bases from the projectors (2.17)–(2.18), self-test |
| `spectral_tools.py` | the diagnostics: unfolding, gap ratio, Wigner surmises, SFF (2.25), cross SFF (2.26), C_AB (2.27) |
| `run_paper_numerics.py` | computes every number and figure at the manuscript's settings |
| `export_couplings.py` | writes the couplings of an ensemble to JSON, for the Mathematica code |
| `crosscheck_export.py` | exports an ensemble for the cross-check against the Mathematica code |
| `S3_rep_category_L8_spectra.npz` | the stored June 2026 ensemble (17 MB, optional) |
| `results/` | output of the November ensemble (the default run) |
| `results_n120_seed20248601/` | output of the regenerated June ensemble |
| `results_S3_rep_category_L8_spectra/` | output of the run on the stored June file |

---

## The three ensembles

|  | Seeds | Reproduce with | Used in the manuscript for |
|---|---|---|---|
| **June 2026** | numpy, base 20248601 | `--seed-base 20248601`, or the stored June file | Table 1, Figures 3 and 4 |
| **November 2026** | numpy, base 20261103 (default) | the default run | Figures 5 and 6 as they now appear |
| **Mathematica** | `SeedRandom[1]` … `SeedRandom[120]` | `../S3_mathematica_paper/` | an independent statistical check |

Realization r draws its couplings with `numpy.random.default_rng(seed_base + r)`, for r = 0, …, 119. The same
couplings enter all three twists, so sectors are paired by realization, as eqs. (2.26) and (2.27) require.

**June.** The June run recorded the seed rule 20240601 + 1000 L + r, which is 20248601 + r at L = 8. This code
draws the couplings in the same order (J, f, J2, J3) with the same means and widths, and builds the same
Hamiltonian and sector bases. This code can therefore make every plot from the June ensemble by regenerating it
at seed base 20248601. The regenerated spectra match the stored June file to 4.5 × 10⁻¹³ in every level of all
120 realizations. The June file itself stores only spectra, not seeds or couplings.

**November.** On the machine below, the default run reproduces the manuscript's Figures 5 and 6 byte for byte.

| File | SHA-256 |
|---|---|
| `results/fig5_sff_GUE.png` | `7b4acab60f09f71d9f5b7b90aea2eb84550634e50f2a29edb38486ed9d0ab3f9` |
| `results/fig6_correlation_matrix.png` | `b1f6f6ea05dae4a2d70fc15c1cab07ce90d1f26a639af7545ac2bf766f918f12` |

> **Caveat.** Byte identity needs the same numpy, LAPACK and matplotlib builds. Elsewhere, expect eigenvalues to
> agree to about 10⁻¹² and the figures to agree visually. Regenerating an ensemble from its seeds relies on
> numpy's random stream, which numpy does not formally guarantee across versions. The stored ensemble files are
> therefore the permanent record.

---

## What reproduces what

| Manuscript | Code | Output |
|---|---|---|
| block dimensions, 1093 to 3281 (3⁸ = 6561 per twisted Hilbert space) | `sector_basis`, `sector_dimension_formula` | `run_log.txt` |
| Table 1, ⟨r⟩ per sector from raw gaps | `rbar_ensemble` | `table1_gap_ratios.json`, `.png` |
| eq. (2.22), Spec H_([c],ω) = Spec H_([c],ω²) | spectra comparison | `run_log.txt` |
| Figure 3, P(s) per sector | `unfold`, `surmise` | `fig3_spacing.png` |
| Figure 4, connected SFF of the six GOE sectors | `amplitudes`, `sff_connected`, `sff_rmt` | `fig4_sff_GOE.png` |
| Figure 5, SFF and cross SFF of the Θ pair | `cross_sff` | `fig5_sff_GUE.png` |
| Figure 6, cross-sector correlation matrix C_AB | `corr_matrix` | `fig6_correlation_matrix.png` |
| "form factors follow RMT for t ≳ 0.2 t_H" | band averages of K_c/K_RMT | `run_log.txt`, `numerics_summary.json` |

---

## Results

All runs of 2026-11-10. Every ensemble has the block dimensions quoted in the manuscript, listed in the table
below. Each twisted Hilbert space has 6561 = 3⁸ states, with ([e],2) counted twice.

**⟨r⟩ per sector (Table 1).** Mean over 120 realizations, with the standard error in brackets.

| Sector | dim | Manuscript | June | November | Mathematica | Class |
|---|---:|---:|---:|---:|---:|:---:|
| ([e],1) | 1094 | 0.5257(10) | 0.5257(10) | 0.5273(11) | 0.5240(11) | GOE |
| ([e],sgn) | 1093 | 0.5262(10) | 0.5262(10) | 0.5261(10) | 0.5257(11) | GOE |
| ([e],2) | 2187 | 0.5289(8) | 0.5289(8) | 0.5274(7) | 0.5291(8) | GOE |
| ([t],+) | 3281 | 0.5298(6) | 0.5298(6) | 0.5303(6) | 0.5294(6) | GOE |
| ([t],−) | 3280 | 0.5288(6) | 0.5288(6) | 0.5298(5) | 0.5303(6) | GOE |
| ([c],1) | 2187 | 0.5294(8) | 0.5294(8) | 0.5292(8) | 0.5286(8) | GOE |
| ([c],ω) | 2187 | 0.5956(7) | 0.5956(7) | 0.5954(7) | 0.5954(8) | GUE |
| ([c],ω²) | 2187 | 0.5956(7) | 0.5956(7) | 0.5954(7) | 0.5954(8) | GUE |

The June column reproduces the manuscript in every digit. The November and Mathematica ensembles agree with it
within 1.4 and 1.7 combined standard errors. The large-matrix benchmarks are 0.5307 (GOE) and 0.5996 (GUE).

**Other diagnostics.** Both codes give the same numbers at the precision shown.

|  | June | November | Mathematica |
|---|:---:|:---:|:---:|
| Θ pair, largest \|E(ω) − E(ω²)\| | 3 × 10⁻¹³ | 3 × 10⁻¹³ | 3 × 10⁻¹³ |
| K_c/K_RMT averaged over 0.2 ≤ t ≤ 1 | 0.95–1.04 | 0.98–1.02 | 0.97–1.02 |
| K_c/K_RMT averaged over 0.02 ≤ t ≤ 0.05, largest → smallest block | 2.1 → 3.8 | 2.0 → 3.6 | 2.1 → 3.9 |
| C_AB of the Θ pair | 1.0000 | 1.0000 | 1.0000 |
| max \|C_AB\| over the 27 distinct pairs (rms) | 0.024 (0.013) | 0.029 (0.014) | 0.032 (0.013) |

The excess at early times is the non-universal regime described in the text. The off-diagonal C_AB are of the
size of the sampling noise of 120 realizations.

---

## Validation

| Check | Result |
|---|---|
| self-test, `python3 s3_chain.py` | all checks pass |
| regenerated June ensemble vs the stored June file | 4.5 × 10⁻¹³ in every level, Table 1 exact |
| per-realization r of the June file vs this code | 5 × 10⁻⁷, their six-decimal rounding |
| default run vs the manuscript's Figures 5 and 6 | byte-identical images (hashes above) |
| Mathematica Hamiltonian on the same couplings | spectra agree to 3.4 × 10⁻¹³ |
| Mathematica diagnostics on identical spectra | ⟨r⟩ 1.5 × 10⁻¹⁵, unfolding 1.6 × 10⁻¹³, SFF 1.7 × 10⁻¹⁰, C_AB 2.5 × 10⁻¹² |
| independent Mathematica ensemble vs November | within 2.2 combined standard errors |

To repeat the Mathematica cross-check, with `OUTDIR` any scratch folder (16 MB):

```bash
python3 crosscheck_export.py OUTDIR
wolframscript -file ../S3_mathematica_paper/crosscheck_with_python.wls OUTDIR
```

---

## Usage in detail

### Commands

```bash
# reproduce the manuscript
python3 run_paper_numerics.py --seed-base 20248601                        # -> results_n120_seed20248601/
python3 run_paper_numerics.py                                             # -> results/

# redraw from a stored ensemble, in seconds
python3 run_paper_numerics.py --from-npz S3_rep_category_L8_spectra.npz   # -> results_S3_rep_category_L8_spectra/
python3 run_paper_numerics.py --from-npz results/ensemble_L8_n120.npz     # -> results/

# explore
python3 run_paper_numerics.py --nreal 10                                  # -> results_n10/
python3 run_paper_numerics.py --seed-base S                               # -> results_n120_seedS/
python3 run_paper_numerics.py ... --out DIR                               # any of the above, -> DIR/
```

* Times were measured on an Apple M4 (macOS 15.6) with Python 3.9.6, numpy 2.0.2, scipy 1.13.1 and
  matplotlib 3.9.4. A full run takes 13 to 19 minutes, depending on the load.
* `--L`, `--fraction`, `--degree`, `--times` and `--band` change the chain length and the analysis settings.
  A value out of range stops the script with an error. `DOCUMENTATION.md` explains them, and how to change the
  disorder strengths.
* Runs that compute another ensemble write to their own folder, so they never overwrite `results/`. A redraw
  with `--from-npz` and default settings writes next to the stored file.
* `--from-npz` reads two formats. One is the format written by this script. The other is the June format, with
  one array per sector and realization, named like `e_1__17`.

### Output files

| File | Contents |
|---|---|
| `run_log.txt` | everything the run printed (`redraw_log.txt` for `--from-npz`) |
| `table1_gap_ratios.json`, `.png` | Table 1 and its plot |
| `fig3_spacing.png` … `fig6_correlation_matrix.png` | Figures 3 to 6 |
| `numerics_summary.json` | the main results and the settings used |
| `ensemble_L8_n120.npz` | all sector spectra, couplings, seed base, numpy version and platform (16 MB) |

### Settings

The manuscript states L = 8, the representative twists e, σ, r, 120 realizations, the central 70% of each
spectrum, degree-8 unfolding, and ⟨r⟩ from the raw gaps. It does not state the time grid, and says only "a band
of unfolded times". This code evaluates the form factors on 80 log-spaced times from 0.02 to 3 Heisenberg times,
and averages C_AB over 0.1 ≤ t ≤ 1. These are the settings of Figures 5 and 6 as they now appear in the
manuscript. Figure 4 used a finer, binned grid (see the notes below).

### Self-test

`python3 s3_chain.py` checks:

* the group action, and U(g) H_h U(g)⁻¹ = H_{ghg⁻¹}
* that H_h commutes with exactly the centralizer of h, and the trace identity
* the orthonormality and invariance of every sector basis, and the dimensions against the character formula
* the recombination of the sector spectra into each full spectrum, and the isospectrality of conjugate twists
* two independent constructions of the ([e],2) block, and the action of Θ on the ([c],ω) basis

It uses `assert`, so do not run it with `python -O`.

---

## Notes on the figures

* **Colours.** The figures follow the manuscript's matplotlib scheme: tab10 colours C0 to C7 for the sectors in
  the order of Table 1, black surmises (Poisson dotted, GOE dashed, GUE dash-dotted), the GOE ramp in black,
  the GUE ramp in red, and the coolwarm map for C_AB.
* **Figure 3 normalization.** The histograms are normalized on 0 ≤ s ≤ 3.5, as in the manuscript. Between
  0.02% and 0.2% of the unfolded gaps lie beyond 3.5, depending on the sector. Many sit in the lowest tenth of
  the window, where the polynomial unfolding is least accurate.
* **Figure 3 legends.** They count all gaps in the central window, for example 91,800 for ([e],1). The
  manuscript shows a few fewer (91,798), because the June script counted only gaps with s < 5 after its own
  unfolding. With this code's unfolding the same cut gives counts within 8 of the June ones. No quoted number
  depends on these counts.
* **Figure 4 time grid.** The manuscript's Figure 4 used 400 log-spaced times from 10^−1.7 to 10^0.5, averaged
  in 44 logarithmic bins. Its curves are therefore smoother than those drawn here from 80 unbinned times.

---

## Packaging for a supplement

The `.npz` ensemble files (about 16 MB each) can be left out, because the scripts regenerate them. Keep the June
file if you want a record of the June ensemble that does not rely on numpy's random stream. Leave out
`__pycache__/`, which Python creates on the first run. On macOS, create the archive with `zip -X` or with
`COPYFILE_DISABLE=1 tar`, so that Finder metadata is not included.
