# $S_3$ qutrit chain

`python/` holds the code behind the $S_3$ numerics of Sec. 2.4.3 (Table 1 and Figs. 3 to 6), written by Sridip Pal. It is included as he wrote it.

| File | Role |
|---|---|
| `python/s3_chain.py` | the model: the $S_3$ action, the twisted Hamiltonian `eq:S3ham`, the disorder `eq:S3disorder`, and the bases of the eight $D(S_3)$ sectors. Run on its own, it is a self-test. |
| `python/spectral_tools.py` | the diagnostics: unfolding, gap ratio, Wigner surmises, and the form factors and correlation matrix of `eq:sff`, `eq:crosssff` and `eq:corrmatrix`. |
| `python/run_paper_numerics.py` | computes Table 1 and Figs. 3 to 6, either from a stored ensemble (`--from-npz`, seconds) or from fresh disorder draws (about 10 minutes). |
| `python/README.md`, `python/DOCUMENTATION.md` (also `.tex` and `.pdf`) | the author's own guide to the package, with every option. |

## How this repository differs from the author's guide

The guide describes the package as it was first assembled. In this repository:

- The two stored ensembles are in `data/s3/`. The guide's `S3_rep_category_L8_spectra.npz` (the June 2026 ensemble) is `data/s3/S3_rep_category_L8_spectra.npz`, and its `results/ensemble_L8_n120.npz` (the November 2026 ensemble) is `data/s3/ensemble_L8_n120.npz`. `data/README.md` describes both.
- The guide's `results*/` output folders are not included. The scripts in `code/figure-reproduction/` regenerate their contents in `code/figure-reproduction/generated/`.
- The Mathematica version of the code (`../S3_mathematica_paper/` in the guide) is not included, and neither are `export_couplings.py` and `crosscheck_export.py`, which prepare the cross-check against it.
- The timings in the guide are for an Apple M4 with Python 3.9.6 and numpy 2.0.2. This repository was tested with the environment of `environment/README.md`, Python 3.14.3 and numpy 2.4.6, where a fresh ensemble takes about 9 minutes.

## Commands

From the repository root, with the interpreter of `environment/README.md`:

```sh
.venv/bin/python code/s3/python/s3_chain.py
.venv/bin/python code/s3/python/run_paper_numerics.py --from-npz data/s3/S3_rep_category_L8_spectra.npz --out code/figure-reproduction/generated/s3_june
.venv/bin/python code/s3/python/run_paper_numerics.py --from-npz data/s3/ensemble_L8_n120.npz --out code/figure-reproduction/generated/s3_nov
```

The first is the self-test. The second gives Table 1 and Figs. 3 and 4 from the June ensemble, and the third gives Figs. 5 and 6 from the November ensemble. Without `--out`, `run_paper_numerics.py` writes to a `results*/` folder inside `code/s3/python/`, which git ignores. The wrappers `s3_from_stored_ensembles.py`, `fig04_s3_sff_goe_binned.py` and `check_s3_regeneration.py` in `code/figure-reproduction/` run these commands and compare the outputs with the paper.
