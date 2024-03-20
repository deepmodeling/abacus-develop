# Solving the Hamiltonian

## Explicit Diagonalization

Method of explicit solving KS-equation can be chosen by variable "ks_solver" in INPUT file.

When "basis_type = pw", `ks_solver` can be `cg`, `bpcg` or `dav`. The `bpcg` method only supports K-point parallelism currently. The default setting `cg` is recommended, which is band-by-band conjugate gradient diagonalization method. There is a large probability that the use of setting of `dav` , which is block Davidson diagonalization method, can be tried to improve performance.  

When "basis_type = lcao", `ks_solver` can be `genelpa` or `scalapack_gvx`. The default setting `genelpa` is recommended, which is based on ELPA (EIGENVALUE SOLVERS FOR PETAFLOP APPLICATIONS) (https://elpa.mpcdf.mpg.de/) and the kernel is auto choosed by GENELPA(https://github.com/pplab/GenELPA), usually faster than the setting of "scalapack_gvx", which is based on ScaLAPACK(Scalable Linear Algebra PACKage)  

## Stochasic DFT
We support stochastic DFT calculation (SDFT) or mixed stochastic-deterministic DFT (MDFT) with plane-wave basis [[Phys. Rev. B 106, 125132 (2022)](https://doi.org/10.1103/PhysRevB.106.125132)]. Different from traditional KSDFT with the explicit diagonalization method, SDFT and MDFT calculate physical quantities with trace of the corresponding operators. The advantages of SDFT and MDFT compared to the traditional KSDFT are the ability to simulate larger sizes and higher temperatures. In our package, SDFT and MDFT can be used by setting the `esolver_type` parameter to `sdft` for SCF calculations or MD calculations. To start with, you can refer to two [examples](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/stochastic) and an explanation of the [input variables](../input_files/input-main.md#electronic-structure-sdft).

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
