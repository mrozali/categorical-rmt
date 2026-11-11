# Analytic checks

This file is secondary context. The paper in `paper/categorical_rmt.tex` is the ground truth.

The checks below were made on 2026-11-10 by Claude (Anthropic), working for the authors while this repository was prepared. They re-derive or test the analytic statements that carry the paper's main results. They are a consistency check, not a peer review. A statement that is not listed here was not checked. Section numbers and equation labels refer to `paper/categorical_rmt.tex`.

## Lattice models (Section 2)

- **Block form.** The block decomposition `eq:Hblocks` follows from Schur's lemma on the irreducible tube-algebra modules $V_\mu$, as stated.
- **$D(S_3)$ sectors (Sec. 2.4).** The eight sectors $([h],R)$ have $\dim V_\mu = |[h]|\dim R = 1, 1, 2, 3, 3, 2, 2, 2$, and the squares add up to $36 = |S_3|^2$. All eight are self-dual, as the footnote in Sec. 2.4.2 says. For $([c],\omega)$ the dual is $(r^2, \bar\omega)$, which conjugation by $\sigma$ carries back to $([c],\omega)$.
- **Multiplicities at $L=8$.** Counting with characters, a transposition fixes one of the $3^8$ configurations and a three-cycle fixes none. This gives multiplicity dimensions $1094$, $1093$ and $2187$ for $([e],1)$, $([e],\mathrm{sgn})$ and $([e],2)$. It gives $3281$ and $3280$ for $([t],\pm)$, and $2187$ for each of the three $[c]$ sectors. That is the range $1093$ to $3281$ quoted in Sec. 2.4.3.
- **$S_3$ classes.** In `eq:S3crt` to `eq:S3degen`, $\Theta$ is complex conjugation in the qutrit basis. It fixes six sectors and exchanges $([c],\omega)$ with $([c],\omega^2)$. So there are six GOE sectors and one GUE pair with coincident spectra.
- **Fibonacci center (Sec. 2.5).** The center $Z(\mathrm{Fib}) = \mathrm{Fib}\boxtimes\overline{\mathrm{Fib}}$ is consistent with the table of $\dim W^\mu_a$ and with `eq:fibdecomp`. The spins of $(\tau,\mathbf 1)$ and $(\mathbf 1,\tau)$ are $e^{\pm 4\pi i/5}$. The actions of $\Theta$ and $\mathcal R\Theta$ in `eq:fibpair` and `eq:fibcrt` also agree with it.
- **Fibonacci block sizes.** The block sizes in `data/fib/fib_L26_lam1.25_spectra.npz` match the dimension statements of Sec. 2.5.2 at $L=26$.
  - The twisted $(\tau,\tau)$ blocks add up to $F_{27} = 196418$ over all 26 momenta.
  - $(\tau,\mathbf 1)$ and $(\mathbf 1,\tau)$ each add up to $F_{26} = 121393$.
  - The file keeps 12 of the 26 untwisted momenta, one from each $\pm k$ pair, with the self-conjugate momenta dropped. Those $W_+$ blocks have 2871 to 2900 levels, against $F_{25}/26 \approx 2886$. The $W_-$ blocks have 7549 to 7560 levels, against $F_{27}/26 \approx 7555$.

## Continuous symmetries (Section 4)

- **$u(1)_k$ conventions.** With $Q = 2J_0$ and $[J_m,J_n] = \tfrac k2 m\,\delta_{m+n,0}$, the Cartan normalization is $K_{U(1)} = 2k$ and $L^{\rm Sug}_n = \tfrac1k\sum :J_mJ_{n-m}:$. The flow formulas (`eq:u1-spectral-flow-algebra_1`, `eq:u1-spectral-flow-algebra_2`, `spectral-flow-state`) agree with the general flow `eq:general-flow`. The vertex operator `equ:SF_vertex_opt` has charge $2k\zeta$ and weight $k\zeta^2$.
- **$u(1)_k$ characters.** The extended characters `2k-theta_function` are index-$k$ theta functions. On them we checked the elliptic law `U(1)_elliptic`, the modular law `U(1)_modular`, the flowed character `eq:flowed-character-general`, and the $S$-dual form `U(1)_modular_eq2`. In `U(1)_modular`, the right-moving sign comes from $\tau_R = -\bar\tau$. The Gaussian $\mathcal G_{U(1)}$ cancels the modular prefactor exactly.
- **$su(2)_k$.** The flow `su(2)_J3_SF` and `su(2)_Jpm_SF` includes a Weyl reflection. It preserves the current algebra and gives $\sigma_\zeta(L_n) = L_n - \zeta J^3_n + \tfrac{k\zeta^2}{4}\delta_{n,0}$. The cocharacter lattice is $2\mathbb Z$. The leading terms of the Weyl–Kac formula `su(2)-char` reproduce $\sum_m e^{2\pi i z m}$ and $h_j - c^{\rm Sug}/24$.
- **Classes (Tables 2 and 3).** $\mathsf{CRT}_V^2 = (-1)^k$ on $V^{\rm ext}_{[k]}$ for $u(1)_k$ and $\mathsf{CRT}_V^2 = (-1)^{2j}$ for $su(2)_k$, as derived in Secs. 4.4.2 and 4.5. Combined with `equ:CRTV-CRTM`, these give both tables.
- **$su(2)_1$ against $u(1)_1$.** The identity $\chi^{SU(2)_1}_j(\tau,z) = \chi^{U(1)_1}_{2j}(\tau,z/2)$ holds, with $Q = 2J^3_0$ and weight $(J^3_0)^2$. The ramp weights $2D^2 = 2, 8, 8, 32$ of the four sectors agree with the $u(1)_1$ count, where the two ground states of charge $\pm1$ give $d_1 = 2$.

## Holography (Sections 3 and 5), read for consistency

- **Spins of $D(S_3)$.** `eq:spins` gives spin $1/3$ to $([c],\omega)$. In a CFT, $\mathsf{CRT}$ fixes this self-dual sector, so it is GOE at fixed spin. Parity pairs it with $([c],\omega^2)$ at opposite spin. This is consistent with Sec. 2, where the disordered chain keeps only $\Theta$ and the pair is GUE.
- **Seed.** The identity kernel is invariant under `eq:Kconj`, because $\varrho(\tilde\gamma) = \overline{\varrho(\gamma)}$. The Fourier transform of the seed reproduces `equ:RMT_result_cate`. It uses $\int dx\, e^{-2\pi i s x}/(x^2+Y^2) = (\pi/Y)\,e^{-2\pi|s|Y}$ with $Y = y_1 + y_2$ and the prefactor $e^{\pi(y_1+y_2)/6}/\sqrt{y_1y_2}$.
- **Parity involution.** The involution $\sigma$ of `eq:sigma` fixes every sector of spin $0$ or $1/2$ in both examples and exchanges the others in pairs, as stated.
- **Sec. 5.**
  - The $c=1$ Kloosterman sums `eq:U1-c1-correction` and `eq:SU2-c1-correction` follow from the stated $S$-matrices.
  - The $SU(2)$ ramp coefficient $D^2/\pi$ in `eq:SU(2)_RMT_result` matches the boundary prediction $2D^2$ times the GUE value $1/(2\pi)$.
  - The GSE slope $\tfrac{1}{4\pi}\times 2^2 = \tfrac1\pi$ equals the GOE slope, as stated.
  - The $SU(2)$ elliptic law `eq:SU2-elliptic` is that of index $k/4$.
- **Fig. 9 (`sf_fig`).** The plotted coordinates agree with `eq:U1-paired-threshold` for $k=1$. The kink sits at $(3/4,\ 13\pi/6)$ for $(Q_L,Q_R) = (2,1)$. The script `code/figure-reproduction/fig09_u1_threshold_check.py` compares every vertex.
- **Appendices.** The flowed weight in App. A.3 equals the weight of spin $k/2 - j$. The identity $J^-_2\,[J^-_{-1}|j,j\rangle]^1 = (k+2j)\,[|j,j\rangle]^1$ holds. The Poincaré-series covariance `eq:PS_VMF` and the Bruhat unfolding behind `eq:PS_Fourier` are standard and consistent with `eq:preK`.

## Not checked

- The full unfolding behind the pre-Kuznetsov formulas `eq:preK` and `WH_amp_result_prime`, beyond the points above.
- The regularization of the gravitational modular sum, which the paper takes from Cotler and Jensen.
- Statements quoted from the literature, such as the Chern–Simons wavefunctions `eq:CS-wavefunctions`.
- The estimates in the Discussion.
