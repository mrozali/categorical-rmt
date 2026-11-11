# Code

- `figure-reproduction/` holds one script for each reproducible figure, plus quick checks of the paper's numerical claims. Its `README.md` maps every figure and table of the paper to a command and a status.
- `fib/` holds the exact-diagonalization code for the perturbed Fibonacci golden chain of Sec. 2.5.
- `s3/` holds the code for the disordered $S_3$ qutrit chain of Sec. 2.4.3. Its `README.md` describes it.

## `fib/`

| File | Role |
|---|---|
| `fib_eth_lowmem.py` | Untwisted chain. Builds the Hamiltonian $H_0 + H_1$ of Sec. 2.5.2 on the fusion-path basis and resolves it by lattice momentum and by the eigenvalue of $\mathcal L^{\tau,\tau}_{\mathbf 1\leftarrow\mathbf 1}$ ($W_\pm$). Diagonalizes each block. Usage: `fib_eth_lowmem.py L lambda [k ...]`. |
| `fib_eth_twisted_lowmem.py` | $\tau$-twisted chain. Resolves it by the characters of the $\mathbb Z_{5L}$ twisted translation, whose $L$-th power gives the topological spin and hence the center label, and diagonalizes each block. Usage: `fib_eth_twisted_lowmem.py L lambda [block ...]`. |
| `fib_chain_L22.py` | numba kernels for the sparse Hamiltonian and for the rows of the $\tau$-line operator, used by both programs. |
| `fib_syk_lib.py` | the chain basis, the $F$-symbols and the twisted-chain helpers, used by the twisted program. |
| `eth_lib.py` | local operators and their matrix elements in the energy basis. |
| `extract_fib_L26_spectra.py` | rebuilds `data/fib/fib_L26_lam1.25_spectra.npz` from the outputs of the two programs at $L=26$. Usage: `extract_fib_L26_spectra.py <output directory>`. |

The two programs also compute matrix elements of local operators, for the companion work on categorical eigenstate thermalization cited in the paper. Only the energy eigenvalues are used here. Both programs write their `.npz` output to the current working directory. Some docstrings and comments mention the authors' internal notes (`anyonic_syk_notes.tex`, `fib_eth_note.tex`) and earlier exploration scripts (paths beginning `fib_chain/`). None of these are part of this repository, and none is needed to run the code.

Blocks with more than 46341 levels would need a 64-bit-integer LAPACK. That code path is not included, and no block used here comes close to it, the largest having 7561 levels.
