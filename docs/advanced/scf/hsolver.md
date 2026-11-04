# Solving the Hamiltonian

## Explicit Diagonalization

Method of explicit solving KS-equation can be chosen by variable "ks_solver" in INPUT file.

When "basis_type = pw", `ks_solver` can be `cg`, `bpcg`, `dav`, `dav_subspace`, or `ppcg`. The default setting `cg` is recommended, which is a band-by-band conjugate-gradient diagonalization method. The `dav` and `dav_subspace` settings use Davidson-style subspace diagonalization and can be tried to improve performance.

The `ppcg` setting uses the projection preconditioned conjugate-gradient method, a restarted block method with bounded workspace. It is optimized and validated for CPU plane-wave calculations; non-CPU devices use a transitional host/device bridge. PPCG is unavailable for `esolver_type = sdft`, including mixed stochastic-deterministic calculations.

PPCG uses `pw_diag_ndim` as its block size and performs the conservative full Rayleigh--Ritz correction on every iteration. Cached Hamiltonian and overlap images are refreshed independently, with additional refreshes when numerical cancellation threatens the residual.

If PPCG reaches `pw_diag_nmax` before satisfying its residual criterion, the running log reports the active-band count and maximum residual. An SCF calculation can continue to its next electronic iteration with a warning; an NSCF calculation stops with an error because it has no subsequent SCF iteration to improve these eigenpairs. A non-finite residual is an error in either calculation.

For diagnostics, `ABACUS_PPCG_RESIDUAL_TRACE` selects an optional CSV trace. Only the MPI world root writes this file, and its resolved full path is recorded in the running log each time the file is created or overwritten.

When "basis_type = lcao", `ks_solver` can be `genelpa` or `scalapack_gvx`. The default setting `genelpa` is recommended, which is based on ELPA (EIGENVALUE SOLVERS FOR PETAFLOP APPLICATIONS) (https://elpa.mpcdf.mpg.de/) and the kernel is auto choosed by GENELPA(https://github.com/pplab/GenELPA), usually faster than the setting of "scalapack_gvx", which is based on ScaLAPACK(Scalable Linear Algebra PACKage)  

## Stochasic DFT
We support stochastic DFT calculation (SDFT) or mixed stochastic-deterministic DFT (MDFT) with plane-wave basis [[Phys. Rev. B 106, 125132 (2022)](https://doi.org/10.1103/PhysRevB.106.125132)]. Different from traditional KSDFT with the explicit diagonalization method, SDFT and MDFT calculate physical quantities with trace of the corresponding operators. The advantages of SDFT and MDFT compared to the traditional KSDFT are the ability to simulate larger sizes and higher temperatures. In our package, SDFT and MDFT can be used by setting the `esolver_type` parameter to `sdft` for SCF calculations or MD calculations. To start with, you can refer to two [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/23_sdft) and an explanation of the [input variables](../input_files/input-main.md#electronic-structure-sdft).

When we have a hamiltonian, the electronic density can be calculated with:

$\rho(\mathbf{r})={\rm Tr}[f(\hat{H})\ket{\mathbf{r}}\bra{\mathbf{r}}]$,

where the Fermi-Dirac function $f(\hat{H})=\frac{1}{1+\exp(\frac{\hat{H}-\mu}{kT})}$ and it can be calculated with the Chebyshev expansion. Here we only support the "fd" or "fermi-dirac" `smearing_method`, the parameter `smearing_sigma` is equal the temperature $T$ (in Ry) and `nche_sto` represents the order of the expansion.

For physical quantities represented by operator $\hat{O}$, SDFT calculates its trace with:

${\rm Tr}[\hat{O}]=\sum_{i=1}^{N_\chi}{\bra{\chi_i}\hat{O}\ket{\chi_i}}$,

while MDFT calculates the trace as:

${\rm Tr}[\hat{O}]=\sum_{n=1}^{N_\phi}{\bra{\phi_n}\hat{O}\ket{\phi_n}}+\sum_{i=1}^{N_\chi}{\bra{\tilde \chi_i}\hat{O}\ket{\tilde \chi_i}}$,

where $\{\ket{\tilde\chi_i}\}$ are obtaiend by projecting stochastic orbitals onto the subspace orthogonal to KS orbitals $\{\phi_n\}$:

$\ket{\tilde\chi_i}=\ket{\chi_i}-\sum_{n=1}^{N_\phi}\braket{\phi_n|\chi_i}\ket{\phi_n}$.

Here the number of KS orbitals $N_\phi$ is controlled by the parameter `nbands` while the number of stochastic orbitals $N_\chi$ is controlled by `nbands_sto`.

Besides, although SDFT does not diagonalize the hamiltonian, it can also caluclate DOS and electronic conductivities with parameters `out_dos` and `cal_cond` separately.
