# Categorical RMT

Chuanxin Cui (University of British Columbia), Kristan Jensen (University of Victoria), Sridip Pal (Institut des Hautes Etudes Scientifiques and Tata Institute of Fundamental Research) and Moshe Rozali (University of British Columbia)

Random matrix theory describes the level statistics of a chaotic system one symmetry sector at a time. This paper proposes its analogue for local quantum systems on a circle with ordinary or non-invertible symmetries, where the sectors are the irreducible representations of the tube algebra. It tests the proposal in chaotic Fibonacci anyon and $S_3$ chains, and extends it to two-dimensional CFTs and their holographic duals.

arXiv: [2611.13093](https://arxiv.org/abs/2611.13093). The paper is also in [`paper/`](paper/), with the compiled PDF at [`paper/categorical_rmt.pdf`](paper/categorical_rmt.pdf). That is the authors' revised version for the second arXiv version, with minor corrections to v1 that [`paper/README.md`](paper/README.md) lists.

## Talk to this paper

This paper is published with an AI agent, following the [Agentic Publication Protocol](https://github.com/LionSR/AgenticPublicationProtocol). Clone this repository and open it in an AI coding agent to ask questions about the paper, check its numbers and reproduce its figures.

**Claude Code:** clone and open the repository, and it reads `AGENTS.md` automatically. Or use `/load-paper https://github.com/mrozali/categorical-rmt`.

**Codex or other agents:** clone and open the repository. Any agent that reads `AGENTS.md` picks up the paper context.

## What is here

- `paper/` holds the LaTeX source, bibliography, figure files and compiled PDF. This is the ground truth.
- `code/fib/` holds the exact-diagonalization code for the Fibonacci chain, and `code/s3/` the code for the $S_3$ chain. `code/README.md` describes both.
- `code/figure-reproduction/` holds the figure scripts and quick checks of the paper's numerical claims.
- `data/` holds the $L=26$ Fibonacci spectra behind Fig. 7 and the two $S_3$ ensembles behind Table 1 and Figs. 3 to 6, about 36 MB in all.
- `environment/` gives the tested platform and the pinned Python packages.
- `supplementary/` holds a record of which analytic statements were re-derived, the validation report and a test of the paper's agent. It is secondary to the paper.

## Figures

The authoritative map from each figure and table to its script and status is [`code/figure-reproduction/README.md`](code/figure-reproduction/README.md).

## Reproducing results

### Setup

[`environment/README.md`](environment/README.md) gives the tested platform, the dependency file, the setup commands and the command prefix for running scripts. In short, from the repository root:

```sh
python3.14 -m venv .venv
.venv/bin/pip install -r environment/requirements.txt
.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py
```

The last command redraws Fig. 7 from the cached spectra in a few seconds and prints the gap ratios quoted in Sec. 2.5.2. `./verify.sh` runs all such checks in one go, in about half a minute, and GitHub Actions runs it on every push.

### Run figures

Use the commands in [`code/figure-reproduction/README.md`](code/figure-reproduction/README.md). Every script writes to `code/figure-reproduction/generated/`.

### Full computation

The cached Fibonacci spectra were computed on a cluster, and recomputing them at $L=26$ needs cluster-scale memory. The programs are in `code/fib/`, and the section "Regenerating the $L=26$ spectra" of the figure map gives the commands. The same programs run at $L=14$ in seconds through `code/figure-reproduction/check_fib_generators_smallL.py`.

Each $S_3$ ensemble can be recomputed from its seeds on a laptop in about 10 minutes, through `code/figure-reproduction/check_s3_regeneration.py`.

## Data

[`data/README.md`](data/README.md) describes the cached spectra, their block labels and how they were made. No download is needed.

## Citation

Please cite the paper. `CITATION.cff` holds the same information in machine-readable form.

```bibtex
@article{Cui:2026catrmt,
    author = "Cui, Chuanxin and Jensen, Kristan and Pal, Sridip and Rozali, Moshe",
    title = "{Categorical RMT}",
    eprint = "2611.13093",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    year = "2026"
}
```

## License

The whole repository (paper, code, data and supplementary material) is released under the Creative Commons Attribution 4.0 International license. See [`LICENSE`](LICENSE).
