# Figure and table reproduction

This file maps every figure and table of `paper/categorical_rmt.tex` to the command that regenerates or checks it. Statuses use the terms of the Agentic Publication Protocol.

**Conventions.** Commands run from the repository root with the interpreter set up in `environment/README.md`, for example `.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py`. Each script writes only to `code/figure-reproduction/generated/`, which is not committed. Each script prints a short summary that ends with `RESULT: PASS` or `RESULT: FAIL`. `./verify.sh` at the repository root runs every script that works from the stored data, in about half a minute.

`s3_from_stored_ensembles.py` is a grouped wrapper. It runs `code/s3/python/run_paper_numerics.py` on both stored $S_3$ ensembles, in about 7 s, and covers Table 1 and Figs. 3, 5 and 6. Its outputs are `generated/s3_june/` and `generated/s3_nov/`, each holding `table1_gap_ratios.json` and `.png`, `fig3_spacing.png`, `fig4_sff_GOE.png`, `fig5_sff_GUE.png`, `fig6_correlation_matrix.png`, `numerics_summary.json` and `redraw_log.txt`. The paper uses Table 1 and Fig. 3 from `s3_june/`, and Figs. 5 and 6 from `s3_nov/`.

## Figures and tables

| Figure/Table | Paper artifact | Script | Inputs | Generated output | Status | Notes |
|---|---|---|---|---|---|---|
| Fig. 1 (`fig:tube`) | TikZ in `paper/categorical_rmt.tex` | none | none | compiled paper | `reproduced` | Schematic, drawn by compiling the paper. |
| Fig. 2 (`fig:group`) | `paper/fig2_group.pdf` | none | none | none | `manual-only` | Drawn schematic, no source file. |
| Table 1 (`tab:S3r`) | table in `paper/categorical_rmt.tex` | `s3_from_stored_ensembles.py` (grouped) | `data/s3/S3_rep_category_L8_spectra.npz` | `generated/s3_june/table1_gap_ratios.json`, `.png` | `reproduced` | June ensemble. Every $\langle r\rangle$ and standard error equals the paper's to the four printed decimals, compared with the values parsed from the paper source. |
| Fig. 3 (`fig:spacing`) | `paper/S3_rep_category_L8_spacing_hist.png` | `s3_from_stored_ensembles.py` (grouped) | `data/s3/S3_rep_category_L8_spectra.npz` | `generated/s3_june/fig3_spacing.png` | `reproduced` | June ensemble. Visually identical. The image is 10 pixels wider, so it is compared by eye. The legend gap counts differ by up to 8 out of about $10^5$ (91,800 against 91,798 for $([e],1)$), because the authors' June plotting script, not included, counted only gaps below $s=5$ after its own unfolding. |
| Fig. 4 (`fig:sffGOE`) | `paper/S3_rep_category_L8_sff_GOE.png` | `fig04_s3_sff_goe_binned.py` | `data/s3/S3_rep_category_L8_spectra.npz` | `generated/fig04_s3_sff_GOE_binned.png` | `runs-but-differs` | June ensemble. The paper's figure was drawn by the authors' June plotting script, which is not included. The script redraws it on the grid that script used, 400 times in 44 logarithmic bins, with a reconstructed binning convention. Same shape and levels, but not point for point. `run_paper_numerics.py` also writes an unbinned 80-time version, `generated/s3_june/fig4_sff_GOE.png`. |
| Fig. 5 (`fig:sffGUE`) | `paper/S3_rep_category_L8_sff_GUE.png` | `s3_from_stored_ensembles.py` (grouped) | `data/s3/ensemble_L8_n120.npz` | `generated/s3_nov/fig5_sff_GUE.png` | `reproduced` | November ensemble. Differs from the paper's image only by anti-aliasing: a largest channel difference of 0.008 on 0.05% of the pixels. Byte-identical on the machine that made it. |
| Fig. 6 (`fig:corr`) | `paper/S3_rep_category_L8_cross_sector.png` | `s3_from_stored_ensembles.py` (grouped) | `data/s3/ensemble_L8_n120.npz` | `generated/s3_nov/fig6_correlation_matrix.png` | `reproduced` | November ensemble. Pixel-identical to the paper's image. |
| Fig. 7 (`numerics_fib`) | `paper/fib_sectors_L26.pdf` | `fig07_fib_sectors_L26.py` | `data/fib/fib_L26_lam1.25_spectra.npz` | `generated/fib_sectors_L26.pdf`, `.png` | `reproduced` | Identical to the paper's figure in a 50 dpi render (0 differing pixels). Also prints the five $\langle r\rangle$ values of Sec. 2.5.2. The $L=26$ spectra are cached, and computing them afresh is `blocked-heavy-compute` (see below). |
| Fig. 8a (`fig:sandwich`) | `paper/fig8a.pdf`, labels in TikZ | none | none | compiled paper | `manual-only` | Drawn schematic, no source file. Its labels are TikZ in the paper source. |
| Fig. 8b, 8c (`fig:sandwich`) | TikZ in `paper/categorical_rmt.tex` | none | none | compiled paper | `reproduced` | Schematics, drawn by compiling the paper. |
| Fig. 9 (`sf_fig`) | TikZ in `paper/categorical_rmt.tex` | `fig09_u1_threshold_check.py` | `paper/categorical_rmt.tex` | terminal | `reproduced` | Recomputes the six plotted vertices and the kink label from `eq:U1-paired-threshold` and compares them with the TikZ coordinates. Max deviation $3.3\times10^{-7}$, the rounding of the coordinates. |
| Table 2 (`table:U1_classes`) | table in `paper/categorical_rmt.tex` | none | none | none | `manual-only` | Analytic, derived in Sec. 4.4.2. Checked by hand, see `supplementary/analytic-checks.md`. |
| Table 3 (`table:SU2_classes`) | table in `paper/categorical_rmt.tex` | none | none | none | `manual-only` | Analytic, derived in Sec. 4.5. Checked by hand, see `supplementary/analytic-checks.md`. |
| Unnumbered tables (Secs. 2.3, 2.4.1, 2.5, 4.3) | tables in `paper/categorical_rmt.tex` | none | none | none | `manual-only` | Dictionaries and sector lists, analytic. The $D(S_3)$ and Fibonacci entries are checked in `supplementary/analytic-checks.md`. |

## Quick claim checks

| Claim | Paper | Command | Expected observation | Evidence level |
|---|---|---|---|---|
| $\langle r\rangle = 0.5292$, $0.5309$, $0.5309$, $0.5306$, $0.5306$ for $W_+$, $W_-$, $(\tau,\tau)$, $(\tau,\mathbf 1)$, $(\mathbf 1,\tau)$ at $L=26$ | Sec. 2.5.2 | `.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py` (about 3 s) | `W+: 0.529223 +- 0.0014` and so on, each line ending `ok` | reproduced from cached spectra |
| $W_-$ and the twisted $(\tau,\tau)$ block have the same spectrum level by level, and so do $(\mathbf 1,\tau)$ and $(\tau,\mathbf 1)$ at opposite momenta, "to within $7\times10^{-14}$" | Sec. 2.5.2 | `.venv/bin/python code/figure-reproduction/check_fib_isospectrality.py` (under 1 s) | `max\|dE\| = 1.066e-14` over 12 pairs for $W_-$ against $(\tau,\tau)$. `max\|dE\| = 6.484e-14 at pair m=97/33` over 26 pairs for $(\tau,\mathbf 1)$ against $(\mathbf 1,\tau)$ | reproduced from cached spectra. Both are coincidences at machine precision. arXiv v1 stated $2\times10^{-14}$, which the largest difference, in one of the 26 pairs, exceeds by about a factor of three. The bound was corrected for the version in `paper/`. |
| The two generators produce the spectra used here and the same coincidences | Sec. 2.5.2 | `.venv/bin/python code/figure-reproduction/check_fib_generators_smallL.py` (about 3 s) | at $L=14$: `5.329e-15` over 6 pairs and `6.217e-15` over 14 pairs | newly checked at small $L$. The paper's claim is at $L=26$. |
| Six of the eight $D(S_3)$ sectors are GOE and $([c],\omega)$, $([c],\omega^2)$ are GUE, with the $\langle r\rangle$ values of Table 1 | Sec. 2.4.3, Table 1 | `.venv/bin/python code/figure-reproduction/s3_from_stored_ensembles.py` (about 7 s) | each Table 1 line ending `ok`, for example `([c],omega)    0.5956 +- 0.0007 \| 0.5956 +- 0.0007 GUE ok` | reproduced from the stored June spectra. The class of each sector also follows analytically from the action of $\Theta$ (Sec. 2.4.2), checked in `supplementary/analytic-checks.md`. |
| The gap ratios lie below the RMT benchmarks "by 1.5 to 5.3 standard deviations" | Sec. 2.4.3, footnote to Table 1 | same command | `1.5 to 5.3 from the unrounded values, 1.5 to 5.7 from the rounded Table 1 entries` | reproduced from the unrounded values. arXiv v1 stated 1.5 to 5.7, which comes from the rounded entries of $([c],\omega)$, $(0.5996-0.5956)/0.0007$. |
| The spectra of the $\Theta$ pair coincide, their cross form factor tracks the auto form factor, $C_{AB}=1$ for the pair and close to 0 for the 27 other pairs | Sec. 2.4.3, Figs. 5 and 6 | same command | June: `Theta pair max \|dE\| = 3.30e-13`, `C_AB(Theta pair) = 1.0000`, `other 27 pairs max \|C_AB\| = 0.024`. November (Fig. 6): `2.83e-13`, `1.0000`, `0.029` | reproduced from the stored spectra, both ensembles |
| The form factors follow RMT for $t\gtrsim0.2\,t_H$, and lie above the ramp at earlier times, most for the smallest blocks | Sec. 2.4.3, Figs. 4 and 5 | `.venv/bin/python code/figure-reproduction/fig04_s3_sff_goe_binned.py` (about 7 s) | $K_c/K_{\rm GOE}$ averaged over $0.2\le t\le1$ between 0.984 and 1.016. Over $0.02\le t\le0.05$, from 2.01 (3280 levels) to 3.56 (1094 levels) | reproduced from the stored June spectra. `s3_from_stored_ensembles.py` prints the same band averages for both ensembles and all eight sectors. |
| The stored $S_3$ spectra follow from the stated model and seeds | Sec. 2.4.3, `eq:S3ham`, `eq:S3disorder` | `.venv/bin/python code/figure-reproduction/check_s3_regeneration.py june` (about 9 min, or `nov`) | `largest difference 4.48e-13` for `june`, `4.02e-13` for `nov` | newly checked, fresh diagonalization of all 120 realizations |
| Fig. 9 follows `eq:U1-paired-threshold` | Sec. 5.2 | `.venv/bin/python code/figure-reproduction/fig09_u1_threshold_check.py` (instant) | `max deviation = 3.33e-07` | reproduced |

## Labels in the Fibonacci data

The stored blocks, their keys and the momentum labels are described in `data/README.md`. Fig. 7 uses the blocks as follows, as Sec. 2.5.2 describes.
- $W_\pm$ uses the 12 stored momenta $k = 1,\dots,12$, one from each $\pm k$ pair.
- $(\tau,\tau)$ uses $m = 5, 10, \dots, 60$. These are the twisted partners of $k = 1, \dots, 12$, so $m = 0$ and the self-conjugate $m = 65$ are dropped.
- $(\tau,\mathbf 1)$ and $(\mathbf 1,\tau)$ use all 26 momenta each.

Each block is trimmed by 5% at each edge and unfolded with a degree-7 polynomial, and the spacings and gap ratios are then pooled.

## Regenerating the $L=26$ spectra

The cached file `data/fib/fib_L26_lam1.25_spectra.npz` was made with the two programs in `code/fib/`. Running them at $L=26$ is `blocked-heavy-compute` for a workstation. The untwisted Hilbert space has $271443$ states and the twisted one $439204$, and every momentum block of up to about 7600 levels is diagonalized densely. The original runs were cluster jobs on Fir, of the Digital Research Alliance of Canada, in July 2026. The commands were:

```sh
.venv/bin/python code/fib/fib_eth_lowmem.py 26 1.25
FIB_PER_BLOCK=1 .venv/bin/python code/fib/fib_eth_twisted_lowmem.py 26 1.25 <block index>
```

The second command was run once for each of the 78 twisted blocks, with block index 0 to 77. Each of those jobs had 48 cores, a memory request of 150 GB and a six-hour limit. These are the requested resources. The peak use was not recorded here, and neither were the resources of the single untwisted job. Each program writes its `.npz` output to the current directory. `code/fib/extract_fib_L26_spectra.py <output directory>` then rebuilds the cached file from those outputs.

Smaller sizes fit on a workstation. The same programs run at $L=14$ in seconds through `check_fib_generators_smallL.py`, and on the tested Mac the twisted program took about 200 s at $L=20$. At $L = 14$ to $20$ they also agreed to about $10^{-14}$ with the authors' earlier, higher-memory implementation of the same chain, which is not included here.

## Regenerating the $S_3$ ensembles

The stored $S_3$ spectra in `data/s3/` can be recomputed from their seeds on a laptop. The command

```sh
.venv/bin/python code/figure-reproduction/check_s3_regeneration.py june
```

runs `code/s3/python/run_paper_numerics.py --seed-base 20248601`. That builds the Hamiltonian `eq:S3ham` for the 120 disorder realizations `default_rng(20248601 + r)` of `eq:S3disorder`, diagonalizes all eight sectors in each, writes Table 1 and Figs. 3 to 6 to `generated/s3_regen_june/`, and compares every level with `data/s3/S3_rep_category_L8_spectra.npz`. The argument `nov` does the same for the November ensemble, with seed base 20261103 and `data/s3/ensemble_L8_n120.npz`. Each takes about 9 minutes on the tested Mac.

On 2026-11-11, with the environment of `environment/README.md`, the regenerated spectra matched the stored ones to $4.5\times10^{-13}$ (June) and $4.0\times10^{-13}$ (November) in every level. The regenerated June ensemble gives Table 1 exactly. From the regenerated November ensemble, Fig. 6 is again pixel-identical to the paper's, and Fig. 5 differs only by anti-aliasing. Regeneration relies on numpy's random stream, which numpy does not formally guarantee across versions, so the stored files are the permanent record.
