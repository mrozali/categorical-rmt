# Paper-agent test

This file is secondary context. It records a test of how an AI agent that opens this repository answers a reader's questions.

## How the test was run

On 2026-11-11 a fresh agent (Claude, Sonnet model, as a general-purpose coding agent) was given a clean copy of this repository as its working directory and nothing else: no conversation history, no access to files outside the copy, and no virtual environment. It was told to start from `AGENTS.md`, to set up the environment from the repository's own docs if it needed to run code, and not to start any computation the docs describe as heavy or longer than a couple of minutes. It then answered twelve reader questions. Its answers were checked against the paper by the agent that prepared the repository (Claude, Opus model).

The repository tested was the release candidate of version 1.0.1. It differs from the release only by this file and `validation-report.md`. Earlier tests of the same kind, on earlier candidates, gave the same answers to the questions they shared.

**Setup it performed.** It read `AGENTS.md`, `README.md`, `code/figure-reproduction/README.md`, `environment/README.md`, `paper/README.md`, `verify.sh` and `CITATION.cff`. It then ran `python3.14 -m venv .venv` and `.venv/bin/pip install -r environment/requirements.txt`, which installed the pinned versions.

**Commands it ran,** all from the repository root:
- `.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py`
- `.venv/bin/python code/figure-reproduction/check_fib_isospectrality.py`
- `.venv/bin/python code/figure-reproduction/s3_from_stored_ensembles.py`
- `./verify.sh`

All ended with `RESULT: PASS` or, for `verify.sh`, `0 failed`.

## Questions and answers

The answers are given as the agent gave them, shortened. The key output lines are kept.

**1. What in this repository should I treat as authoritative, and what is only secondary context?**

> The paper, the code and the data are the ground truth: `paper/categorical_rmt.tex` (with its PDF), `code/fib/` and `code/s3/python/`, the figure map `code/figure-reproduction/README.md` and its scripts, and `data/` with `data/README.md`. `environment/README.md` is the source of truth for setup. Secondary: `supplementary/analytic-checks.md`, a consistency check made for the authors, the validation records in `supplementary/`, and the change note `paper/README.md`. Not included: the companion paper on categorical eigenstate thermalization and its data, the authors' June plotting script, and some internal notes named in code comments. Keep three kinds of support apart. The sector structure, the symmetry classes and Tables 2 and 3 are analytic. The lattice results are finite-size numerics. The CFT predictions are not tested directly, and the holography reaches the ramp but not the plateau.

**2. In one paragraph, what is the paper's main claim, and what does it predict that ordinary symmetry-resolved random matrix theory would not?**

> The paper proposes "categorical RMT". For a local chaotic system on a circle with an ordinary or non-invertible symmetry, the independent random-matrix blocks are the sectors of the tube algebra, not the twisted Hilbert spaces. The Hamiltonian on each multiplicity space, restricted to any further unitary symmetry such as lattice momentum, is an ordinary Wigner-Dyson matrix whose class is fixed by the antiunitary symmetries that map the sector to itself, and distinct blocks are independent unless a symmetry relates them (Sec. 2.2). For a group this reproduces symmetry-resolved RMT. For a non-invertible symmetry, Hilbert spaces with different twists can share a sector, and so part of their spectrum, level by level, although no invertible symmetry relates them. In the Fibonacci chain the untwisted W₋ block and the twisted (τ,τ) block are different Hamiltonians on spaces of different dimension, yet have identical levels. The paper also predicts GSE neutral sectors in bosonic CFTs with U(1)_k or SU(2)_k symmetry, and shows that decorated torus wormholes reproduce the ramp and its modular completion, though they cannot tell GOE from GSE.

**3. The paper quotes mean gap ratios ⟨r⟩ for the five Fibonacci blocks at L=26. Can you check them from what is here, and tell me what you get?**

```
.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py
  W+: 0.529223 +- 0.0014  rounds to 0.5292, expected 0.5292  ok
  W-: 0.530853 +- 0.0009  rounds to 0.5309, expected 0.5309  ok
  W1: 0.530853 +- 0.0009  rounds to 0.5309, expected 0.5309  ok
  W2: 0.530572 +- 0.0008  rounds to 0.5306, expected 0.5306  ok
  W3: 0.530572 +- 0.0008  rounds to 0.5306, expected 0.5306  ok
pixel check (50 dpi, ...): max diff 0.0000, fraction of differing pixels 0.00000
RESULT: PASS
```

> Yes, from the cached spectra (not a fresh diagonalization). The five values match Sec. 2.5.2 to the four printed decimals and lie close to ⟨r⟩_GOE = 0.5307. The identical values for W₋ and (τ,τ), and for (τ,1) and (1,τ), reflect the level-by-level coincidences. The redrawn Fig. 7 has 0 differing pixels against the paper's figure.

**4. Sec. 2.5.2 gives a numerical bound on how closely the shared-sector levels coincide. Does the data support it?**

```
.venv/bin/python code/figure-reproduction/check_fib_isospectrality.py
(i)  W-(k) vs (tau,tau)(m=5k), 12 pairs: max|dE| = 1.066e-14   (paper states 7e-14)
(ii) (tau,1)(m) vs (1,tau)(130-m), 26 pairs: max|dE| = 6.484e-14 at pair m=97/33   (paper states 7e-14)
RESULT: PASS
```

> Yes. Both families lie below the stated 7×10⁻¹⁴, the second with about 8% to spare. arXiv v1 stated 2×10⁻¹⁴, which the second family exceeds, and `paper/README.md` records the correction. This is machine-precision agreement of the stored eigenvalues. The L=14 check of the same programs gives about 5×10⁻¹⁵ and 6×10⁻¹⁵ with fresh spectra, but the paper's claim is at L=26.

**5. Can I reproduce Table 1 (the S3 gap ratios) from this repository? If so, do it and tell me how it compares.**

```
.venv/bin/python code/figure-reproduction/s3_from_stored_ensembles.py
  ([e],1)        0.5257 +- 0.0010 | 0.5257 +- 0.0010 GOE ok | 0.5273 +- 0.0011
  ...
  ([c],omega^2)  0.5956 +- 0.0007 | 0.5956 +- 0.0007 GUE ok | 0.5954 +- 0.0007
  November against June: largest difference 1.5 combined standard errors
RESULT: PASS
```

> Yes, from the stored June ensemble in about 7 s. Every ⟨r⟩ and standard error equals the paper's to the four printed decimals, six sectors GOE and the pair ([c],ω), ([c],ω²) GUE. The same run reproduces Figs. 3, 5 and 6 (Fig. 6 pixel-identical, Fig. 5 differing only by anti-aliasing). I did not run the 9-minute regeneration from the seeds, which the repository records as passing.

**6. Are Table 1 and Figs. 3 to 6 all computed from the same disorder realizations?**

> No, and Sec. 2.4.3 says so: there are two independent sets of 120 realizations at L=8. The June ensemble (`data/s3/S3_rep_category_L8_spectra.npz`, seeds `default_rng(20248601 + r)`) gives Table 1 and Figs. 3 and 4. The November ensemble (`data/s3/ensemble_L8_n120.npz`, seeds `20261103 + r`) gives Figs. 5 and 6. Their gap ratios agree within 1.5 combined standard errors, and in both the off-diagonal C_AB stay below 0.03. arXiv v1 described a single set. The June file records no seeds, and its seed rule comes from an earlier script that is not included.

**7. Is the footnote's statement about how far the S3 gap ratios lie below the RMT benchmarks right?**

```
offset below the RMT benchmark, in standard errors: 1.5 to 5.3 from the unrounded values, 1.5 to 5.7 from the rounded Table 1 entries   (paper footnote: 1.5 to 5.3, arXiv v1: 1.5 to 5.7)
```

> Yes, in this version: the unrounded values give exactly 1.5 to 5.3. arXiv v1 said 5.7, from the rounded ([c],ω) entry. The footnote's explanation, weak disorder breaking reflection only slightly, is the authors' suggestion and is not tested here.

**8. Why can a neutral sector have symplectic (GSE) statistics in a bosonic CFT with U(1)_k symmetry? Point me to where the paper shows this.**

> On a self-conjugate sector CRT factorizes as CRT_V ⊗ CRT_V̄ ⊗ CRT_M, and the class is set by CRT_M². In a bosonic theory CRT² = +1, so CRT_M² = ε_L ε_R with ε_L = CRT_V² and ε_R = CRT_V̄² (`equ:CRTV-CRTM`). In Sec. 4.4.2 the charge-k module contains |k⟩ and |−k⟩, related by the spectral-flow operator of spin k. CRT maps a bosonic operator of spin s to (−1)^s times its conjugate (`equ:CRT_SF_opt`), which gives CRT_V² = (−1)^k. For (k,0) or (0,k) with k odd the product is −1, so the neutral block is GSE with Kramers pairs (Table 2). The SU(2)_k analogue, half-integer j_L+j_R, is Table 3. The torus wormhole gives these sectors the GOE slope, equal to the GSE slope with Kramers pairs counted twice, so holography cannot tell them apart. This is an analytic derivation, not tested numerically.

**9. I want to recompute the L=26 Fibonacci spectra from scratch on my laptop. What will that take?**

> It is not feasible on a laptop, and the repository marks it `blocked-heavy-compute`. The untwisted space has 271,443 states and the twisted one 439,204, with dense diagonalization of blocks up to about 7,600 levels. The original runs were one untwisted job and 78 twisted jobs on the Fir cluster, each twisted job requesting 48 cores, 150 GB and six hours. What fits on a laptop: the same programs at L=14 in seconds (`check_fib_generators_smallL.py`), and the twisted program at L=20 in about 200 s per the docs.

**10. I have the arXiv v1 PDF. Is the paper in this repository the same, and if not, what changed?**

> Not identical. Per `paper/README.md`, this is the authors' revised version for the second arXiv version, with the same section numbers and these changes:
> 1. the reference `Pal:2020wwd`, which showed as `[?]` in v1,
> 2. the bound in Sec. 2.5.2, now 7×10⁻¹⁴ instead of 2×10⁻¹⁴,
> 3. a new sentence in Sec. 2.4.3 on the two sets of realizations,
> 4. the footnote range, now 1.5 to 5.3 instead of 5.7,
> 5. typos in Secs. 2.5, 4.1, 4.3, 4.4.2, 5.1 and 6.
>
> I confirmed items 2, 3 and 4 from the source and the scripts.

**11. Is there a single command that checks everything that can be checked quickly? Run it and tell me the result.**

```
./verify.sh
verify.sh quick, interpreter .venv/bin/python
OK    s3_selftest (0 s)
OK    fig07_fib_sectors_L26 (1 s)
OK    check_fib_isospectrality (0 s)
OK    check_fib_generators_smallL (5 s)
OK    fig09_u1_threshold_check (0 s)
OK    s3_from_stored_ensembles (6 s)
OK    fig04_s3_sff_goe_binned (5 s)
0 failed
```

> Yes, `./verify.sh` from the repository root: 0 failures in about 15 s. It covers Fig. 7 and the five ⟨r⟩ values, the isospectrality bounds, the L=14 generator check, Fig. 9, Table 1, and Figs. 3 to 6. It does not cover the L=26 recomputation, the fresh S3 diagonalization (`./verify.sh full` adds it), the analytic results or the schematics. Fig. 4 passes its own checks but is `runs-but-differs` against the paper's image.

**12. How should I cite this work?**

> Cite the paper, with the BibTeX entry given in `AGENTS.md` and `README.md` (key `Cui:2026catrmt`, arXiv 2611.13093, hep-th). `CITATION.cff` carries the same paper as its preferred citation and describes the repository itself, version 1.0.1, under CC BY 4.0.

## Assessment

| Question | Tests | Result |
|---|---|---|
| 1 | ground truth and secondary material | correct |
| 2 | main claim, and what is new beyond ordinary RMT | correct |
| 3 | headline Fibonacci numbers from the cached data | ran the check, reported the observed values |
| 4 | the isospectrality bound | ran the check, reported the observed values and the margin |
| 5 | Table 1 | ran the check, every row equal, separated stored spectra from a fresh computation |
| 6 | provenance of the S3 figures | correct, both ensembles identified |
| 7 | the footnote's range | correct. In its full answer it mixed units when quoting the largest offset, though the footnote's comparison it was checking holds. |
| 8 | analytic derivation, with paper anchors | correct, labels verified. It called Sec. 4.3 "Sec. 4.4" once. |
| 9 | heavy computation | warned, did not run it, offered the small-size runs |
| 10 | relation to arXiv v1 | correct |
| 11 | the one-command check | ran it, 0 failed, and stated what it does not cover |
| 12 | citation | correct |

The agent answered from the files in the repository, used commands that run from its root, reported observed values for every numerical question, and stated each limitation the docs record.
