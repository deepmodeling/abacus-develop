# Solving the Hamiltonian

## Explicit Diagonalization

### Plane Wave
#### CG
#### Davidson

### LCAO
#### blablabla

## Stochasic DFT
We support stochastic DFT calculation (SDFT) or mixed stochastic-deterministic DFT (MDFT) with plane-wave basis [[Phys. Rev. B 106, 125132 (2022)](https://doi.org/10.1103/PhysRevB.106.125132)]. Different from traditional KSDFT with the explicit diagonalization method, SDFT and MDFT calculate physical quantities with trace of the corresponding operators. The advantages of SDFT and MDFT compared to the traditional KSDFT are the ability to simulate larger sizes and higher temperatures. In our package, SDFT and MDFT can be used by setting the `calculation` parameter to `sto-scf` or `sto-md` for SCF calculations or MD calculations. To start with, you can refer to an easy [example](../../examples/stochastic.md) and an explanation of the [input variables](../../input-main.md#electronic-structure-sdft).

When we have a hamiltonian, the electronic density can be calculated with 
$$\rho(\mathbf{r})={\rm Tr}[f(\hat{H})\ket{\mathbf{r}}\bra{\mathbf{r}}],$$
where the Fermi-Dirac function $f(\hat{H})=\frac{1}{1+\exp(\frac{\hat{H}-\mu}{kT})}$ and it can be calculated with the Chebyshev expansion. Here we only support the "fd" or "fermi-dirac" `smearing_method`, the parameter `smearing_sigma` is equal the temperature $T$ (in Ry) and `nche_sto` represents the order of the expansion.
SDFT calculates the trace with 
$${\rm Tr}[\hat{O}]=\sum_{i=1}^{N_\chi}{\bra{\chi_i}\hat{O}\ket{\chi_i}},$$
while MDFT calculates the trace with
$${\rm Tr}[\hat{O}]=\sum_{n=1}^{N_\phi}{\bra{\phi_n}\hat{O}\ket{\phi_n}}+\sum_{i=1}^{N_\chi}{\bra{\tilde \chi_i}\hat{O}\ket{\tilde \chi_i}},$$
where $\ket{\tilde\chi_i}$s are orthogonal to KS orbitals $\{\phi_n\}$:
$$\ket{\tilde\chi_i}=\ket{\chi_i}-\sum_{n=1}^{N_\phi}\braket{\phi_n|\chi_i}\ket{\phi_n}.$$
Here the number of KS orbitals $N_\phi$ is controlled by the parameter `nbands` while the number of stochastic orbitals $N_\chi$ is controlled by `nbands_sto`.

Besides, although SDFT does not diagonalize the hamiltonian, it can also caluclate DOS and electronic conductivities with parameters `out_dos` and `cal_cond` separately.
