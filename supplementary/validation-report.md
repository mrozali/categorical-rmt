# Validation report

This file is secondary context. It records the full validation of this repository against the Agentic Publication Protocol (APP) 1.1.0 before its first release.

- **Stage:** full.
- **Result:** passed.
- **Date:** 2026-11-11.
- **Validated by:** Claude (Anthropic, Opus model), working for the authors. The paper-agent test used a separate fresh agent (Sonnet model).
- **Release:** `mrozali/categorical-rmt`, tag `v1.0.0`, version `1.0.0`.
- **Paper:** arXiv 2611.13093. The source in `paper/` is the authors' revised version for the second arXiv version (see `paper/README.md`).

## Checks passed

**Structure.**
- The top level holds only `paper/`, `code/`, `data/`, `environment/`, `supplementary/`, `AGENTS.md`, `CLAUDE.md` (`@AGENTS.md`), `README.md`, `LICENSE` and `.gitignore`.
- There is one canonical paper source, `paper/categorical_rmt.tex`, declared as `latex` in the `AGENTS.md` frontmatter.
- No hidden, cached or generated files are committed. Script outputs, the virtual environment and caches are gitignored.

**Paths and commands.**
- Every path named in `AGENTS.md`, `README.md` and the READMEs under `paper/`, `code/`, `data/` and `environment/` resolves from the repository root.
- Every shell block in those files parses.
- The two external links (the protocol repository and the arXiv abstract) return HTTP 200. The repository URL itself resolves once the repository is public.

**Figure and table map.**
- `code/figure-reproduction/README.md` lists every figure and table of the paper with a final status:
  - `reproduced`: Fig. 1, Table 1, Figs. 3, 5, 6, 7, 8b, 8c and 9.
  - `runs-but-differs`: Fig. 4. The paper's version was drawn by a plotting script that is not included.
  - `manual-only`: the drawn schematics, Fig. 2 and Fig. 8a, and the analytic tables.
- Each non-reproduced item gives its reason.
- Computing the $L=26$ Fibonacci spectra afresh is marked `blocked-heavy-compute`, with the cluster resources the original jobs requested.

**Commands run.** All were run from the repository root with the pinned environment of `environment/README.md` (Python 3.14.3 on macOS, Apple silicon). Every one ended with `RESULT: PASS` (or, for the self-test, `self-test passed`):

| Command | Time |
|---|---|
| `fig07_fib_sectors_L26.py` | about 3 s |
| `check_fib_isospectrality.py` | under 1 s |
| `check_fib_generators_smallL.py` | about 3 s |
| `fig09_u1_threshold_check.py` | instant |
| `s3_from_stored_ensembles.py` | about 7 s |
| `fig04_s3_sff_goe_binned.py` | about 7 s |
| `check_s3_regeneration.py june` | 513 s |
| `code/s3/python/s3_chain.py` (self-test) | under a second |

The scripts without a path are in `code/figure-reproduction/`. The November $S_3$ ensemble was also regenerated with the same code and matched its stored file to $4.0\times10^{-13}$.

**Paper build.** `paper/categorical_rmt.tex` compiles with the four-pass build of `environment/README.md` to 77 pages, with no errors and no undefined references or citations.

**Quick claim checks.** Every numerical statement of Secs. 2.4.3 and 2.5.2 has a quick check in the figure map, giving a command, the expected output and an evidence level. This covers:
- the five Fibonacci gap ratios
- the isospectrality bound
- Table 1
- the footnote to Table 1
- the coincidence of the $\Theta$ pair and $C_{AB}$
- the time range of RMT behaviour in the form factors
- the regeneration of the stored $S_3$ spectra.

**Data.**
- `data/README.md` describes both data sets, their formats, provenance and uses.
- All data are in the repository, about 36 MB, and no download is needed.

**Consistency.**
- Title, authors, arXiv identifier and BibTeX entry agree between `AGENTS.md` and `README.md`.
- `version` 1.0.0 matches the tag `v1.0.0`.
- The paper summary and key results in `AGENTS.md` were checked against the paper source.
- Statements about reproduction status agree across `README.md`, `AGENTS.md`, the figure map and this report.

**Privacy.**
- No credentials, private paths, hostnames or IP addresses appear in the tracked files.
- The only e-mail addresses are the authors' addresses on the paper's title page and one in the header of the third-party `JHEP.bst`.

**License.** `LICENSE` is the CC BY 4.0 legal code, verbatim, and covers the paper, code, data and supplementary material.

**Paper-agent test.** A fresh agent with only this repository answered ten reader questions, recorded in `paper-agent-test.md`. It set up the environment from the docs, ran four checks, and reported observed values for every numerical question. It also stated each limitation that the docs record. All ten answers were correct.

## Limits recorded, not blockers

- Fig. 4 is `runs-but-differs` because the plotting script behind the paper's image is not available. The same spectra and settings give the same curves within noise.
- The $L=26$ Fibonacci spectra are cached. Recomputing them needs cluster-scale memory, while the same programs run at $L=14$ in seconds.
- Figs. 2 and 8a are drawn schematics with no source file.
- `supplementary/analytic-checks.md` lists the analytic statements that were re-derived, and those that were not.

## Public-release blockers

None for APP compliance. Release still requires the authors' approval, a public repository, the tag and the release manifest.
