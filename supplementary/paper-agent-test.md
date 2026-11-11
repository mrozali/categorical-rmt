# Paper-agent test

This file is secondary context. It records a test of how an AI agent that opens this repository answers a reader's questions.

## How the test was run

On 2026-11-11 a fresh agent (Claude, Sonnet model, as a general-purpose coding agent) was given a clean copy of this repository as its working directory and nothing else: no conversation history, no access to files outside the copy, and no virtual environment. It was told to start from `AGENTS.md`, to set up the environment from the repository's own docs if it needed to run code, and not to start any computation the docs describe as heavy or longer than a couple of minutes. It then answered ten reader questions. Its answers were checked against the paper by the agent that prepared the repository (Claude, Opus model). Every value and line reference it quoted was confirmed.

The repository tested was the release candidate, which differs from the published version only by this file and `validation-report.md`.

**Setup it performed.** It read `AGENTS.md`, `code/figure-reproduction/README.md`, `environment/README.md`, `data/README.md` and `paper/README.md`. It then ran `python3.14 -m venv .venv` and `.venv/bin/pip install -r environment/requirements.txt`, which succeeded with the pinned versions.

**Commands it ran,** all from the repository root, all ending `RESULT: PASS`:
- `fig07_fib_sectors_L26.py`
- `check_fib_isospectrality.py`
- `check_fib_generators_smallL.py`
- `s3_from_stored_ensembles.py`

All four are in `code/figure-reproduction/` and were run with `.venv/bin/python`.

## Questions and answers

The answers are given as the agent gave them, lightly shortened. The commands and their key output lines are kept.

**1. What in this repository should I treat as authoritative, and what is only secondary context?**

> The paper, the code and the data are the ground truth. Secondary context is everything else.
>
> Authoritative:
> - `paper/categorical_rmt.tex` is the canonical source, and `paper/categorical_rmt.pdf` is compiled from it. Its equation and section labels (for example `eq:catrmtintro`, `sec:fibnum`) are the stable anchors.
> - `code/fib/` and `code/s3/python/` hold the exact-diagonalization programs.
> - `code/figure-reproduction/README.md` is the authoritative map from every figure and table to its command and status.
> - `data/` holds the cached spectra, and `data/README.md` gives their labels and provenance.
> - `environment/README.md` is the source of truth for setup.
>
> Secondary:
> - `supplementary/analytic-checks.md` is a consistency check made for the authors, not part of the paper.
> - `paper/README.md` records the differences from arXiv v1.
>
> Not included: the companion paper on categorical eigenstate thermalization (`Jensen:2026cat`) and its data, the authors' June plotting script, and some internal notes named in code comments.
>
> Keep three kinds of support apart. The sector structure, the symmetry classes and Tables 2 and 3 are analytic. The lattice results are finite-size numerics. The CFT predictions are not tested directly, and the holography reaches the ramp but not the plateau.

**2. In one paragraph, what is the paper's main claim, and what does it predict that ordinary symmetry-resolved random matrix theory would not?**

> The paper proposes "categorical RMT". For a local chaotic system on a circle with ordinary or non-invertible symmetry, the independent random-matrix blocks are the irreducible sectors of the tube algebra, not the twisted Hilbert spaces. The Hamiltonian on each multiplicity space, restricted to any further unitary symmetry such as momentum, is an ordinary Wigner-Dyson matrix whose class is fixed by the antiunitary symmetries that map the sector to itself, and distinct blocks are independent unless a symmetry relates them. For a group this reproduces ordinary symmetry-resolved RMT. The new content is for non-invertible symmetry: Hilbert spaces with different twists can share a sector, and hence part of their spectrum, level by level, even though no invertible symmetry relates them. In the Fibonacci chain the untwisted and τ-twisted chains are different Hamiltonians on spaces of different dimension, yet W₋ and the twisted (τ,τ) have the same levels. Ordinary RMT would treat the two spaces as unrelated. The paper also predicts that in CFTs with U(1)_k or SU(2)_k symmetry the neutral spectrum can be GSE in a bosonic theory (Sec. 4, Tables 2 and 3), and that torus wormholes decorated by the SymTFT or Chern-Simons theory reproduce the ramp, with the GSE sectors showing the GOE slope. Evidence level: analytic for the structure and classes, numerical at finite size for the Fibonacci and S3 lattice checks (S3 is a group, so it is a check, not a new prediction), and the CFT predictions are not directly tested.

**3. The paper quotes mean gap ratios ⟨r⟩ for the five Fibonacci blocks at L=26. Can you check them from what is here, and tell me what you get?**

```
.venv/bin/python code/figure-reproduction/fig07_fib_sectors_L26.py
  W+: 0.529223 +- 0.0014  rounds to 0.5292, expected 0.5292  ok
  W-: 0.530853 +- 0.0009  rounds to 0.5309, expected 0.5309  ok
  W1: 0.530853 +- 0.0009  rounds to 0.5309, expected 0.5309  ok
  W2: 0.530572 +- 0.0008  rounds to 0.5306, expected 0.5306  ok
  W3: 0.530572 +- 0.0008  rounds to 0.5306, expected 0.5306  ok
pixel check (50 dpi ...): max diff 0.0000, fraction of differing pixels 0.00000
RESULT: PASS
```

> Yes, from the cached L=26 spectra (reproduced from cached spectra, not recomputed). I get W₊ = 0.5292 ± 0.0014, W₋ = 0.5309 ± 0.0009, twisted (τ,τ) = 0.5309 ± 0.0009, (τ,1) = 0.5306 ± 0.0008 and (1,τ) = 0.5306 ± 0.0008. These match Sec. 2.5.2 to the printed four decimals, and all lie within about 1 SEM of ⟨r⟩_GOE = 0.5307. The Fig. 7 render is pixel-identical to the paper's figure. W₋ and (τ,τ) coincide exactly because they are the same spectra (see Q4).

**4. Sec. 2.5.2 gives a numerical bound on how closely the shared-sector levels coincide. Does the data support it?**

```
.venv/bin/python code/figure-reproduction/check_fib_isospectrality.py
(i)  W-(k) vs (tau,tau)(m=5k), 12 pairs: max|dE| = 1.066e-14   (paper states 7e-14)
(ii) (tau,1)(m) vs (1,tau)(130-m), 26 pairs: max|dE| = 6.484e-14 at pair m=97/33   (paper states 7e-14)
RESULT: PASS
.venv/bin/python code/figure-reproduction/check_fib_generators_smallL.py
W-(k) vs (tau,tau)(m=5k): 6 pairs, max|dE| = 5.329e-15
(tau,1)(m) vs (1,tau)(70-m): 14 pairs, max|dE| = 6.217e-15, worst pair m=17/53
RESULT: PASS
```

> Yes. The paper says the coincidences hold "to within 7×10⁻¹⁴". The cached L=26 spectra give 1.07×10⁻¹⁴ over the 12 W₋/(τ,τ) pairs and 6.48×10⁻¹⁴ over the 26 (τ,1)/(1,τ) pairs, both below the bound, with little margin for the second. These are coincidences at machine precision. arXiv v1 said 2×10⁻¹⁴, which the data exceed, and the version here corrects it. This is a check on cached spectra, not a fresh L=26 computation. The L=14 run of the same two generator programs also coincides, which supports the generators but is a smaller size than the paper's claim.

**5. Can I reproduce Table 1 (the S3 gap ratios) from this repository? If so, do it and tell me how it compares.**

```
.venv/bin/python code/figure-reproduction/s3_from_stored_ensembles.py
  ([e],1)        0.5257 +- 0.0010 | 0.5257 +- 0.0010 GOE ok | 0.5273 +- 0.0011
  ...
  ([c],omega^2)  0.5956 +- 0.0007 | 0.5956 +- 0.0007 GUE ok | 0.5954 +- 0.0007
  November against June: largest difference 1.5 combined standard errors
RESULT: PASS
```

> Yes, from the stored June spectra (`data/s3/S3_rep_category_L8_spectra.npz`), in seconds. Every ⟨r⟩ and standard error equals the paper's Table 1 to the four printed decimals. Six sectors are GOE (about 0.526 to 0.530) and ([c],ω), ([c],ω²) are GUE at 0.5956. The same run shows the Θ pair's spectra coinciding to 3.3×10⁻¹³, C_AB = 1.0000 for the pair and at most 0.024 for the other 27 pairs, and reproduces Figs. 5 and 6 from the November ensemble. This confirms the table from the stored spectra, not by fresh diagonalization. The docs report a regeneration check (`check_s3_regeneration.py june`, about 9 minutes, matching to 4.5×10⁻¹³) that I did not run because it exceeds my time limit.

**6. Are Table 1 and Figs. 3 to 6 all computed from the same disorder realizations?**

> No. Sec. 2.4.3 says "two independent sets of 120 disorder realizations each. Table 1 and Figs. 3 and 4 use the first set, and Figs. 5 and 6 the second." In the repository the first set is the June ensemble (`data/s3/S3_rep_category_L8_spectra.npz`, seeds `default_rng(20248601 + r)`), and the second is the November ensemble (`data/s3/ensemble_L8_n120.npz`, seeds `20261103 + r`). The script compares them: the largest ⟨r⟩ difference is 1.5 combined standard errors, and C_AB stays below 0.03 in both. arXiv v1 described a single set. The June file records no seeds, and its seed rule comes from the authors' earlier script, which is not included.

**7. Is the footnote's statement about how far the S3 gap ratios lie below the RMT benchmarks right?**

```
offset below the RMT benchmark, in standard errors: 1.5 to 5.3 from the unrounded values, 1.5 to 5.7 from the rounded Table 1 entries   (paper footnote: 1.5 to 5.3, arXiv v1: 1.5 to 5.7)
```

> Yes, for the version in this repository: from the unrounded values the offsets span 1.5 to 5.3 standard errors. The 5.7 in arXiv v1 comes from the rounded entry of ([c],ω), (0.5996 − 0.5956)/0.0007. The footnote's explanation, weak disorder breaking reflection only slightly, is the authors' suggestion ("a plausible source") and is not tested here.

**8. Why can a neutral sector have symplectic (GSE) statistics in a bosonic CFT with U(1)_k symmetry? Point me to where the paper shows this.**

> This is derived analytically, not tested numerically. The class is set by the square of the antiunitary CRT on the neutral multiplicity space, and the theory is bosonic, so CRT² = +1 on the whole Hilbert space. CRT factorizes as CRT_V ⊗ CRT_V̄ ⊗ CRT_M, so CRT_M² = ε_L ε_R with ε_L = CRT_V² and ε_R = CRT_V̄² (eq. `equ:CRTV-CRTM`, Sec. 4.3). If the two current modules have opposite signs, the neutral part squares to −1 and is GSE, with every neutral level an exact Kramers pair. For U(1)_k (Sec. 4.4.2), the ground states |k⟩ and |−k⟩ of the charge-k module are related by the spectral-flow vertex operator, which has spin k, and CRT maps a bosonic operator of spin s to (−1)^s times its Hermitian conjugate (eq. `equ:CRT_SF_opt`). This gives CRT_V² = (−1)^k on V_[k], so the sectors (k,0) and (0,k) mod 2k with k odd are GSE (Table 2, `table:U1_classes`). For SU(2)_k the analogue is the modules with half-integer j_L + j_R (Table 3). On the holographic side (Sec. 5), the torus wormhole gives these sectors the GOE slope 1/π, equal to the GSE slope 1/(4π) × 2² once each Kramers pair is counted twice, so the wormhole cannot tell GOE from GSE. No numerical test exists, since the paper says no chaotic irrational CFT is under enough control.

**9. I want to recompute the L=26 Fibonacci spectra from scratch on my laptop. What will that take?**

> It is not feasible on a laptop, so I did not attempt it. The figure map marks it `blocked-heavy-compute`. The untwisted Hilbert space has 271,443 states and the twisted one 439,204, and every momentum block of up to about 7600 levels is diagonalized densely. The original runs were cluster jobs (Fir, Digital Research Alliance of Canada, July 2026): one untwisted job, `fib_eth_lowmem.py 26 1.25`, and 78 twisted jobs, `FIB_PER_BLOCK=1 fib_eth_twisted_lowmem.py 26 1.25 <block index>`, each requesting 48 cores, 150 GB of memory and a six-hour limit. `extract_fib_L26_spectra.py` then rebuilds the cached file. What fits on a laptop: the same programs at L=14 in seconds (`check_fib_generators_smallL.py`, which I ran, see Q4), and L=20 for the twisted program in about 200 s on the tested Mac, per the docs.

**10. I have the arXiv v1 PDF. Is the paper in this repository the same, and if not, what changed?**

> No. Per `paper/README.md`, this is the authors' revised version, prepared on 2026-11-11 for the second arXiv version, and it differs from v1 only by minor corrections, with the same section numbers. Changes:
> 1. The missing reference `Pal:2020wwd`, which showed as `[?]` in a footnote to Sec. 2.2 of v1.
> 2. Sec. 2.5.2: the coincidence bound is now 7×10⁻¹⁴, where v1 said 2×10⁻¹⁴ (my Q4 run gives 6.484×10⁻¹⁴).
> 3. Sec. 2.4.3: a new sentence says Table 1 and Figs. 3 and 4 use one set of 120 realizations and Figs. 5 and 6 a second (confirmed in the source).
> 4. Table 1 footnote: 1.5 to 5.3 standard deviations, where v1 had 5.7 (my Q7 run confirms both numbers).
> 5. Typos in Secs. 2.5, 4.1, 4.3, 4.4.2, 5.1 and 6.
>
> To verify against your v1 PDF, compare items 2, 3 and 4 first, since they are the substantive ones.

## Assessment

| Question | Tests | Result |
|---|---|---|
| 1 | ground truth and secondary material | correct |
| 2 | main claim, and what is new beyond ordinary RMT | correct |
| 3 | headline Fibonacci numbers from the cached data | ran the check, reported the observed values |
| 4 | the isospectrality bound | ran both checks, reported the observed values and the margin |
| 5 | Table 1 | ran the check, every row equal, separated "from stored spectra" from "fresh computation" |
| 6 | provenance of the S3 figures | correct, both ensembles identified |
| 7 | the footnote's range | ran the check, reported both the corrected and the v1 value |
| 8 | analytic derivation, with paper anchors | correct, labels and sections verified |
| 9 | heavy computation | warned, did not run it, offered the small-size runs |
| 10 | relation to arXiv v1 | correct. One overstatement in the full answer: it said none of the changes alter a claim, while items 2 and 4 correct stated numbers. |

The agent answered from the files in the repository, used commands that run from its root, reported observed values for every numerical question, and stated each limitation the docs record.
