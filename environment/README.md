# Environment

This file is the source of truth for recreating the computational environment of this repository.

## Summary

- Tested platform: macOS 26.5.2 on Apple silicon (arm64), 2026-11-10.
- Primary toolchain: Python 3.14.3. TeX Live 2025 (pdfTeX and BibTeX) to compile the paper.
- Environment required: yes, a Python virtual environment with five pinned packages.
- External or proprietary requirements: none for the default workflow. Computing the $L=26$ spectra afresh needs cluster-scale memory (see `code/figure-reproduction/README.md`).

## Dependency files

- `environment/requirements.txt` pins the Python packages: numpy 2.4.6, scipy 1.17.1, matplotlib 3.10.9, numba 0.65.1 and llvmlite 0.47.0. llvmlite is a dependency of numba. All other imports come from the Python standard library or from the modules in `code/fib/` and `code/s3/python/`.

## Setup from the repository root

```sh
python3.14 -m venv .venv
.venv/bin/pip install -r environment/requirements.txt
```

Run every script from the repository root with the prefix `.venv/bin/python`, for example:

```sh
.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py
```

Newer releases of the packages were not tried.

`./verify.sh` runs, in one command, the self-test of the $S_3$ code and every script of `code/figure-reproduction/` that works from the stored data, and reports each as OK or FAIL (about half a minute). `./verify.sh full` adds the two $S_3$ regenerations, about 9 minutes each. It uses `.venv/bin/python` when that exists, and the variable `PYTHON` names another interpreter. GitHub Actions runs `./verify.sh quick` on every push, on Ubuntu with Python 3.14 and the same pinned packages (`.github/workflows/verify.yml`).

## Other tools

| Software | Version tested | Used for | Required? |
|---|---|---|---|
| pdfTeX, BibTeX | TeX Live 2025 | compiling `paper/categorical_rmt.tex` | only to rebuild the PDF |
| pdftoppm (Poppler) | 26.04.0 | the pixel comparison at the end of `fig07_fib_sectors_L26.py` | optional, the script skips that comparison without it |

To rebuild the paper, run from `paper/`:

```sh
pdflatex categorical_rmt && bibtex categorical_rmt && pdflatex categorical_rmt && pdflatex categorical_rmt
```

The result has 77 pages, with no errors and no undefined references.

## Gitignored local artifacts

- `.venv/` is the virtual environment.
- `code/fib/__pycache__/` holds numba's compiled-kernel cache, which appears after the first run of the Fibonacci programs.
- `code/figure-reproduction/generated/` holds all script outputs.
- `code/s3/python/results*/` holds the output of `run_paper_numerics.py` when it is run without `--out`.
- LaTeX build files, if the paper is compiled inside `paper/`.

## Validation evidence

| Command | Working directory | Result | Notes |
|---|---|---|---|
| `python3.14 -m venv .venv` and `.venv/bin/pip install -r environment/requirements.txt` | repository root | passed | 2026-11-10, Python 3.14.3, pinned versions installed in about 8 s |
| `.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py` | repository root | `RESULT: PASS` | about 3 s |
| `.venv/bin/python code/figure-reproduction/check_fib_isospectrality.py` | repository root | `RESULT: PASS` | under 1 s |
| `.venv/bin/python code/figure-reproduction/check_fib_generators_smallL.py` | repository root | `RESULT: PASS` | about 3 s |
| `.venv/bin/python code/figure-reproduction/fig09_u1_threshold_check.py` | repository root | `RESULT: PASS` | instant |
| `.venv/bin/python code/s3/python/s3_chain.py` | repository root | `self-test passed` | 2026-11-11, under a second |
| `./verify.sh quick` | repository root | `0 failed` | 2026-11-11, 7 checks in about 15 s |
| `.venv/bin/python code/figure-reproduction/s3_from_stored_ensembles.py` | repository root | `RESULT: PASS` | 2026-11-11, about 7 s |
| `.venv/bin/python code/figure-reproduction/fig04_s3_sff_goe_binned.py` | repository root | `RESULT: PASS` | 2026-11-11, about 7 s |
| `.venv/bin/python code/figure-reproduction/check_s3_regeneration.py june` | repository root | `RESULT: PASS` | 2026-11-11, 513 s, largest difference $4.5\times10^{-13}$. The November ensemble, regenerated with the same code, matched to $4.0\times10^{-13}$. |
| the four-pass build of the paper | copy of `paper/` | passed | 77 pages, 0 errors, 0 undefined references |
