# Constructing the Hamiltonian

## Exchange-Correlation Functionals

In our package, the XC functional can be set explicitly using the `dft_functional` keyword in `INPUT` file. If `dft_functional` is not specified, ABACUS will use the xc functional indicated in the pseudopotential file. 

Several common functionals are implemented in ABACUS, such as PZ and PBE. Users can check out this [file](https://github.com/deepmodeling/abacus-develop/blob/develop/source/source_hamilt/module_xc/xc_funcs.h) for a complete list of functionals implemented in ABACUS. Furthermore, if ABACUS is compiled with LIBXC, we also support all the LDA, GGA and meta-GGA functionals provided therein.

Here, we use a simple [example calculation](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/02_scf/02_lcao_Si2) for illustration.

1. **Default setting:**

    In the original `INPUT` file, there is no specification of the `dft_functional` keyword. As a result, we use the default option, that is to use the xc functional in the pseudopotential file, `Si.pz-vbc.UPF`. We can take a look at the first few lines of the `<PP_HEADER>` section from the pseudopotential file:

    ```
    <PP_HEADER>
    0                   Version Number
    Si                   Element
    NC                  Norm - Conserving pseudopotential
        F                  Nonlinear Core Correction
    SLA  PZ   NOGX NOGC   PZ   Exchange-Correlation functional
    ```

    From the line commented 'Exchange-Correlation functional', we see that this pseudopotential is generated using PZ functional. As a result, if we run ABACUS with the original setting, PZ functional will be used.

    > Note : for systems with multiple elements, if no `dft_functional` is specified, users should make sure that all pseudopotentials are using the same functional. Otherwise, the type of xc functional should be specified explicitly.

2. **Using PBE**

    On the other hand, users might also explicitly specify the xc functional through `dft_functional` parameter. For example, to use PBE functional, add the following line to `INPUT` file and rerun the calculation:

    ```
    dft_functional PBE
    ```

3. **More functionals from LIBXC**

    ABACUS has its own implementation of the PBE functional as well as a few others, but our list is far from comprehensive. However, if ABACUS is compiled with LIBXC, we also support all the LDA, GGA and meta-GGA functionals provided therein.

    For this part, users should compile the ABACUS code with LIBXC linked (version 5.1.7 or higher).

    To use SCAN functional, make the following modification to the `INPUT` file:

    ```
    dft_functional SCAN
    ```
    
    Note that in the case of PBE and SCAN, we are using 'short-hand' names to represent the entire functional, which is made up of individual exchange and correlation components. A complete list of 'short-hand' expressions supported by ABACUS can be found in [source code](../../../source/source_hamilt/module_xc/xc_functional.cpp).

    Apart from the 'short-hand' names, ABACUS also allow supplying exchange-correlation functionals as combinations of LIBXC keywords for functional components, joined by plus sign, for example, setting:

    ```
    dft_functional LDA_X_YUKAWA+LDA_C_1D_CSC
    ```
    means we are using the short-range Yukawa attenuated exchange along with the Casula, Sorella & Senatore LDA correlation functional.

    The list of LIBXC keywords can be found on its [website](https://www.tddft.org/programs/libxc/functionals/).

4. **Temperature-dependent functional**

    In ABACUS, we provide temperature-dependent functionals through LIBXC. For such functionals, the keyword `xc_temperature` (unit is Rydberg) is used to specify the temperature, such as the following:

    ```
    dft_functional LDA_XC_CORRKSDT 
    xc_temperature 10
    ```

5. **Hybrid functional**

    ABACUS supports functionals with exact Hartree-Fock exchange in LCAO basis set only. The old INPUT parameter exx_hybrid_type for hybrid functionals has been absorbed into `dft_functional`. Options are `hf` (pure Hartree-Fock), `pbe0`(PBE0), `hse`, and `scan0`(SCAN0) (Note: in order to use HSE or SCAN0 functional, LIBXC is required). Note also that only HSE has been tested while other hybrid functionals have NOT been fully tested yet, and the maximum parallel cpus for running exx is N^4, with N being the number of atoms.

    More information on the hybrid functional can be found from the section [Exact Exchange](../input_files/input-main.md#exact-exchange-common) in the list of input variables for more information.

    An example HSE calculation is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/20_hybrid/02_lcao_Si2). Apart from the input files (`INPUT`, `STRU`, `KPT`), we further provide two files: running_scf.log_ref and log_ref, which contains reference for running_scf.log and standard output from the program, respectively.
    
## DFT+*U*

Conventional functionals, e.g., L(S)DA and GGAs, encounter failures in strongly correlated systems, usually characterized by partially filled *d*/*f* shells. These include transition metals (TM) and their oxides, rare-earth compounds, and actinides, to name a few, where L(S)DA/GGAs typically yield quantitatively or even qualitatively wrong results. To address this failure, an efficient and successful method named DFT+*U*, which inherits the efficiency of L(S)DA/GGA but gains the strength of the Hubbard model in describing the physics of strongly correlatedsystems, has been developed.

Now the DFT+*U* method is accessible in ABACUS. The details of the DFT+*U* method could be found in this [paper](https://doi.org/10.1063/5.0090122). It should be noted that the DFT+*U* works only within the NAO scheme, which means that the value of the keyword `basis_type` must be lcao when DFT+*U* is called. To turn on DFT+*U*, users need to set the value of the `dft_plus_u` keyword in the `INPUT` file to be 1. All relevant parmeters used in DFT+*U* calculations are listed in the [DFT+*U* correction](../input_files/input-main.md#dftu-correction) part of the [list of keywords](../input_files/input-main.md).

Examples of DFT+*U* calculations are provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/19_dftu).

### Full-projector double counting

Here **Full** names the existing radius-adjustable localized projection
(`dft_plus_u=1`, controlled by `onsite_radius`) implemented by `DFTU_onsite`.
It does **not** mean a full anisotropic $U/J$ interaction. The input
`hubbard_u` is the simplified $U_{\rm eff}=U-J$; no separate Hund's $J$ or
Slater-integral tensor is introduced. Retaining `DFTU_onsite` is not the
addition of a new On-site/projector-3 method.

For this projector, `dftu_dc` selects `fll` (default),
`amf`, or `mixed`. Omitting the new parameters or explicitly selecting `fll`
keeps the existing FLL potential and energy path. This extension does not
change Dual (`dft_plus_u=2`) or add another projector.

For fixed-weight interpolation, specify both parameters:

```text
dftu_dc     mixed
dftu_alpha  0.5
```

`dftu_alpha` must be explicitly supplied, finite and in [0,1] for `mixed`;
it is rejected for other modes. Zero gives AMF and one gives FLL. The mixing
is applied to the local interaction energy and its potential, not to two
separate SCF results. This is a fixed input weight, not a self-consistently
determined interpolation parameter.

For each atom and spin channel $\sigma$, let $n^\sigma$ be the existing
projected Hermitian occupation matrix of dimension $M=2l+1$ and
$N_\sigma=\mathrm{Tr}(n^\sigma)$. AMF subtracts the orbital average
$\bar n_\sigma=N_\sigma/M$ **within each spin channel**, not the common
spin average $(N_\uparrow+N_\downarrow)/(2M)$. No additional trace
normalization or rescaling of the existing projected occupation is applied.
Suppressing the atom/spin labels, the isotropic effective-$U$ convention is

$$
E_\alpha=\frac{U_{\rm eff}}{2}
\left[\alpha N+(1-\alpha)\frac{N^2}{M}-\mathrm{Tr}(n^2)\right],
\qquad
V_\alpha=U_{\rm eff}\left[
\left(\frac{\alpha}{2}+(1-\alpha)\frac{N}{M}\right)I-n\right].
$$

In particular, $E_{\rm AMF}=-U_{\rm eff}\mathrm{Tr}[(n-\bar n I)^2]/2$.
Thus a uniform shell in either spin channel has zero AMF correction, even
when the two spin channels have different occupations. Both the functional
and its derivative are interpolated at the **same** occupation and a fixed
input $\alpha$, not at occupations from separate converged calculations.

Here $\delta E=\mathrm{Re}\,\mathrm{Tr}(V\delta n)$. Energies are summed over
atoms and spin channels. For `nspin=1`, the existing Full occupation builder
halves the spin-summed density; the resulting shell energy is doubled, but
the potential is not. The Full operator's band-energy correction is
$E_\alpha-\mathrm{Re}\,\mathrm{Tr}(V_\alpha n)$, evaluated with the **same**
occupation as the potential. This extension retains that existing update
timing; it does not implement a different post-diagonalization energy ledger.
Equivalently, for a spin-summed matrix $D$ in `nspin=1`,
$E_{\rm total}(D)=2E_\alpha(D/2)$ and
$\delta E_{\rm total}=\mathrm{Re}\,\mathrm{Tr}[V_\alpha(D/2)\delta D]$.
The factor two in the energy cancels the one-half from differentiating the
occupation, so doubling the potential would be incorrect. With `nspin=2`,
the two shell energies and band corrections are accumulated once each.

AMF and mixed requests, including both endpoints, currently require explicit
PBE, double-precision CPU LCAO, fixed-structure KS SCF, collinear `nspin=1/2`,
`symmetry=0` or `-1`, one k pool and atomic cold start. The supported solver
choices are `genelpa`, `elpa`, `lapack` and `scalapack_gvx`. Gamma and multi-k
representations use the same local-shell formulas. SOC/spin constraints,
force/stress requests even during fixed-structure SCF, restart/occupation
files, U ramping, Yukawa U, occupation/DMR
mixing, EXX and matrix/Wannier/QO exports are rejected. These new restrictions
do **not** apply to legacy FLL requests. Charge and Mulliken output remain
available.

For context on FLL, AMF and interpolation, see
[Dudarev et al., Phys. Rev. B 57, 1505 (1998)](https://doi.org/10.1103/PhysRevB.57.1505),
[Czyzyk and Sawatzky, Phys. Rev. B 49, 14211 (1994)](https://doi.org/10.1103/PhysRevB.49.14211), and
[Petukhov et al., Phys. Rev. B 67, 153106 (2003)](https://doi.org/10.1103/PhysRevB.67.153106).
The equations above specify the convention implemented here; this is not a
general Slater-integral or self-consistent-alpha implementation.
