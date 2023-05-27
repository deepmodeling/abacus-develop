# Constructing the Hamiltonian

## Exchange-Correlation Functionals

In our package, the XC functional can be set explicitly using the `dft_functional` keyword in `INPUT` file. If `dft_functional` is not specified, ABACUS will use the xc functional indicated in the pseudopotential file. 

Several common functionals are implemented in ABACUS, such as PZ and PBE. Users can check out this [file](../../../source/module_hamilt_general/module_xc/xc_funcs.h) for a complete list of functionals implemented in ABACUS. Furthermore, if ABACUS is compiled with LIBXC, we also support all the LDA, GGA and meta-GGA functionals provided therein.

Here, we use a simple [example calculation](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/scf/lcao_Si2) for illustration.

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
    
    Note that in the case of PBE and SCAN, we are using 'short-hand' names to represent the entire functional, which is made up of individual exchange and correlation components. A complete list of 'short-hand' expressions supported by ABACUS can be found in [source code](../../../source/module_hamilt_general/module_xc/xc_functional.cpp).

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

    More information on the hybrid functional can be found from the section [Exact Exchange](../input_files/input-main.md#exact-exchange) in the list of input variables for more information.

    An example HSE calculation is provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/hse/lcao_Si2). Apart from the input files (`INPUT`, `STRU`, `KPT`), we further provide two files: running_scf.log_ref and log_ref, which contains reference for running_scf.log and standard output from the program, respectively.
    
## DFT+*U*

Conventional functionals, e.g., L(S)DA and GGAs, encounter failures in strongly correlated systems, usually characterized by partially filled *d*/*f* shells. These include transition metals (TM) and their oxides, rare-earth compounds, and actinides, to name a few, where L(S)DA/GGAs typically yield quantitatively or even qualitatively wrong results. To address this failure, an efficient and successful method named DFT+*U*, which inherits the efficiency of L(S)DA/GGA but gains the strength of the Hubbard model in describing the physics of strongly correlatedsystems, has been developed.

Now the DFT+*U* method is accessible in ABACUS. The details of the DFT+*U* method could be found in this [paper](https://doi.org/10.1063/5.0090122). It should be noted that the DFT+*U* works only within the NAO scheme, which means that the value of the keyword `basis_type` must be lcao when DFT+*U* is called. To turn on DFT+*U*, users need to set the value of the `dft_plus_u` keyword in the `INPUT` file to be 1. All relevant parmeters used in DFT+*U* calculations are listed in the [DFT+*U* correction](../input_files/input-main.md#dftu-correction) part of the [list of keywords](../input_files/input-main.md).

Examples of DFT+*U* calculations are provided in this [directory](https://github.com/deepmodeling/abacus-develop/tree/develop/examples/dft_plus_u).
