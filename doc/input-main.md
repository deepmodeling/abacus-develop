# INPUT file
- [Structure of the file](#structure-of-the-file)
- [List of keywords](#list-of-keywords)
    - [System variables](#system-variables)

        [suffix](#suffix) | [ntype](#ntype) | [nbands](#nbands) | [atom_file](#atom-file) | [kpoint_file](#kpoint-file) | [pseudo_dir](#pseudo-dir) | [nbands_istate](#nbands-istate) | [nspin](#nspin) | [calculation](#calculation) | [dft_functional](#dft-functional) | [read_file_dir](#read-file-dir) | [pseudo_type](#pseudo-type)
    - [Plane wave related variables](#plane-wave-related-variables)
    
        [ecutwfc](#ecutwfc) | [ethr](#ethr) | [start_wfc](#start-wfc) | [start_charge](#start-charge)

    - [Electronic structure and geometry relaxation](#electronic-structure-and-geometry-relaxation)
    
        [basis_type](#basis-type) | [ks_solver](#ks-solver) | [smearing](#smearing) | [sigma](#sigma) | [mixing_type](#mixing-type) | [mixing_beta](#mixing-beta) | [mixing_ndim](#mixing-ndim) | [mixing_gg0](#mixing-gg0) | [gamma_only](#gamma-only) | [printe](#printe) | [niter](#niter) | [diago_cg_maxiter](#diago-cg-maxiter) | [diago_david_ndim](#diago-david-ndim) | [dr2](#dr2) | [charge_extrap](#charge-extrap) | [out_charge](#out-charge) | [out_potential](#out-potential) | [out_dm](#out-dm) | [out_wf](#out-wf) | [out_lowf](#out-lowf) | [mulliken](#mulliken) | [out_alllog](#out-alllog) | [force](#force) | [nstep](#nstep) | [force_thr](#force-thr) | [force_thr_ev](#force-thr-ev) | [bfgs_w1](#bfgs-w1) | [bfgs_w2](#bfgs-w2) | [trust_radius_max](#trust-radius-max) | [trust_radius_min](#trust-radius-min) | [trust_radius_ini](#trust-radius-ini) | [stress](#stress) | [stress_thr](#stress-thr) | [press](#press) | [fixed_axes](#fixed-axes) | [move_method](#move-method)

    - [Molecular dynamics](#molecular-dynamics)

        [md_type](#md-type) | [md_rstmd](#md-rstmd) | [md_dt](#md-dt) | [md_t](#md-t) | [md_qmass](#md-qmass) | [md_nresn](#md-nresn) | [md_dumpmdfred](#md-dumpmdfred) | [md_domsd](#md-domsd) | [md_fixtemperature](#md-fixtemperature) | [md_msdstarttime](#md-msdstarttime)
    - [VdW correction](#vdw-correction)

        [vdw_method](#vdw-method) | [vdw_s6](#vdw-s6) | [vdw_s8](#vdw-s8) | [vdw_a1](#vdw-a1) | [vdw_a2](#vdw-a2) | [vdw_d](#vdw-d) | [vdw_abc](#vdw-abc) | [vdw_C6_file](#vdw-C6-file) | [vdw_C6_unit](#vdw-C6-unit) | [vdw_R0_file](#vdw-R0-file) | [vdw_R0_unit](#vdw-R0-unit) | [vdw_model](#vdw-model) | [vdw_radius](#vdw-radius) | [vdw_radius_unit](#vdw-radius-unit) | [vdw_cn_radius](#vdw-cn-radius) | [vdw_cn_radius_unit](#vdw-cn-radius-unit) | [vdw_period](#vdw-period)
    - [Berry phase and wannier90 interface](#berry-phase-and-wannier90-interface)
    
        [berry_phase](#berry-phase) | [gdir](#gdir) | [towannier90](#towannier90) | [nnkpfile](#nnkpfile) | [wannier_spin](#wannier-spin) | [tddft](#tddft)  [vext](#vext) | [vext_dire](#vext-dire) 

    [back to main page](../README.md)

## Structure of the file
Below is an example INPUT file with some of the most important parameters that need to be set:
```
INPUT_PARAMETERS
#Parameters (General)
ntype 1
nbands 4
#Parameters (Accuracy)
ecutwfc 60
```
Parameters list starts with key word `INPUT_PARAMETERS`. Any content before `INPUT_PARAMETERS` will be ignored.

Any line starting with `#` will also be ignored.

Each parameter value is provided by specifying the name of the input variable
and then putting the value after the name, separated by one or more blank characters(space
or tab). The following characters(≤ 150) in the same line will be neglected.

Depending on the input variable, the value may be an integer, a real number or a string. The parameters can be given in any order, but only one parameter should be given per line.

Furthermore, if a given parameter name appeared more than once in the input file, only the last value will be taken.

*Note: if a parameter name is not recognized by the program, the program will stop with an
error message.*

In the above example, the meanings of the parameters are:
- ntype : how many types of elements in the unit cell
- nbands : the number of bands to be calculated
- ecutwfc : the plane-wave energy cutoff for the wave function expansion (UNIT: Rydberg)

## List of keywords
    
### System variables
This part of variables are used to control general system parameters.
- suffix<a id="suffix"></a>
    - *Type*: String
    - *Description*: In each run, ABACUS will generate a subdirectory in the working directory. This subdirectory contains all the information of the run. The subdirectory name has the format: OUT.suffix, where the ‘suffix’ is the name you can pick up for your convenience.
    - *Default*: ABACUS

    [back to top](#input-file)
- ntype<a id="ntype"></a>
    - *Type*: Integer
    - *Description*: Number of different atom species in this calculations. This value must be set. If the number you set is smaller than the atom species in the STRU file, ABACUS only read the ‘wrong number’ of atom information. If the number is larger than the atom species in the STRU file, ABACUS may stop and quit.
    - *Default*: ***No default value***

    [back to top](#input-file)
- nbands<a id="nbands"></a>
    - *Type*: Integer
    - *Description*: Number of bands to calculate. It is recommended you setup this value, especially when you use smearing techniques, more bands should be included.
    - *Default (nspin=1)*: *1.2\*occupied_bands, occupied_bands+10)*
    - *Default (nspin=2)*: *max(1.2\*nelec, nelec+20)*

    [back to top](#input-file)
- atom_file<a id="atom-file"></a>
    - *Type*: String
    - *Description*: This parameter specifies the name of structure file which contains various information about atom species, including pseudopotential files, local orbitals files, cell information, atom positions, and whether atoms should be allowed to move.
    - *Default*: STRU

    [back to top](#input-file)
- kpoint_file<a id="kpoint-file"></a>
    - *Type*: String
    - *Description*: This parameter specifies the name of k-points file. Note that if you use atomic orbitals as basis, and you only use gamma point, you don’t need to have k-point file in your directory, ABACUS will automatically generate ‘KPT’ file. Otherwise, if you use more than one k-point, please do remember the algorithm in ABACUS is different for gamma only and various k-point dependent simulations. So first you should turn off the k-point algorithm by set `gamma_only = 0` in `INPUT` and then you should setup your own k-points file.
    - *Default*: KPT

    [back to top](#input-file)
- pseudo_dir<a id="pseudo-dir"></a>
    - *Type*: String
    - *Description*: This parameter specifies pseudopotential directory.
    - *Default*: ./

    [back to top](#input-file)
- nbands_istate<a id="nbands-istate"></a>
    - *Type*: Integer
    - *Description*: Only used when `calculation = ienvelope` or `calculation = istate`, this variable indicates how many bands around Fermi level you would like to calculate. `ienvelope` means to calculate the envelope functions of wave functions <em>&Psi;<sub>i</sub>=&Sigma;<sub>&mu;</sub>C<sub>i&mu;</sub>&Phi;<sub>&mu;</sub></em>, where <em>&Psi;<sub>i</sub></em> is the ith wave function with the band index <em>i</em> and <em></sub>&Phi;<sub>&mu;</sub></em> is the localized atomic orbital set. `istate` means to calculate the density of each wave function <em>|&Psi;<sub>i</sub>|<sup>2</sup></em>. Specifically, suppose we have highest occupied bands at 100th wave functions. And if you set this variable to 5, it will print five wave functions from 96th to 105th. But before all this can be carried out, the wave functions coefficients  should be first calculated and written into a file by setting the flag `out_lowf = 1`.
    -   *Default*: 5

    [back to top](#input-file)
- nspin<a id="nspin"></a>
    - *Type*: Integer
    - *Description*: Number of spin components of wave functions. There are only two choices now: 1 or 2, meaning non spin or collinear spin.
    - *Default*: 1

    [back to top](#input-file)
- calculation<a id="calculation"></a>
    - *Type*: String
    - *Description*: Specify the type of calculation.
        - *scf*: do self-consistent electronic structure calculation
        - *relax*: do structure relaxation calculation, one can ues ‘nstep’ to decide how many ionic relaxations you want.
        - *cell-relax*: do cell relaxation calculation.
        - *nscf*: do the non self-consistent electronic structure calculations. For this option, you need a charge density file.
        - *istate*: Please see the explanation for variable `nbands_istate`.
        - *ienvelope*: Please see the explanation for variable `nbands_istate`.
        - *md*: molecular dynamics
    - *Default*: scf
    
    [back to top](#input-file)

- dft_functional<a id="dft-functional"></a>
    - *Type*: String
    - *Description*: In our package, the XC functional can either be set explicitly using the dft_functional keyword as explained below, or set implicitly according to the XC functional information read from pseudopotential file. The user should ensure that the XC functional set in the INPUT file and the pseudopotential file are consistent.

        To be specific, we briefly explain the format of the pseudopotential file and the key information it contains. There are a few lines in Si’s GGA pseudopotential file Si_ONCV_PBE-1.0.upf:
        ```
        ...
        <PP_HEADER
        generated="Generated using ONCVPSP code by D. R. Hamann"
        author="Martin Schlipf and Francois Gygi"
        date="150105"
        comment=""
        element="Si"
        pseudo_type="NC"
        relativistic="scalar"
        is_ultrasoft="F"
        is_paw="F"
        is_coulomb="F"
        has_so="F"
        has_wfc="F"
        has_gipaw="F"
        core_correction="F"
        functional="PBE"
        z_valence=" 4.00"
        total_psenergy=" -3.74274958433E+00"
        rho_cutoff=" 6.01000000000E+00"
        ```

        Possible values of this variable are:
        - none: the functional is specified implicity by the input pseudopotential file
        - lda: Perdew-Zunger local density approximation
        - pbe: Perdew-Burke-Ernzerhof general gradient approximation
    
        If the functional specified by the user is not consistent with the pseudopotential file, the program will stop with an error message.
    - *Default*: none

    [back to top](#input-file)

- read_file_dir<a id="read-file-dir"></a>
    - *Type*: String
    - *Description*: when the program needs to read files such as electron density(`SPIN1_CHG`) as a starting point, this variables tells the location of the files. For example, './' means the file is located in the working directory.
    - *Default*: OUT.$suffix

    [back to top](#input-file)

- pseudo_type<a id="pseudo-type"></a>
    - *Type*: String
    - *Description*: the format of pseudopotential files. Accepted value s are:
        - upf : .UPF format
        - vwr : .vwr format
        - upf201 : the new UPF format
    - *Default* : upf

    [back to top](#input-file)

### Plane wave related variables
This part of variables are used to control the plane wave related parameters.

- ecutwfc<a id="ecutwfc"></a>
    - *Type*: Real
    - *Description*: Energy cutoff for plane wave functions, the unit is Rydberg. Note that even for localized orbitals basis, you still need to setup a energy cutoff for this system. Because our local pseudopotential parts and the related force are calculated from plane wave basis set, etc. Also, because our orbitals are generated by matching localized orbitals to a chosen set of wave functions from certain energy cutoff, so this set of localize orbitals are most accurate under this same plane wave energy cutoff.
    - *Default*: 50

    [back to top](#input-file)
- ethr<a id="ethr"></a>
    - *Type*: Real
    - *Description*: Only used when you use diago_type = cg or diago_type = david. It indicates the threshold for the first electronic iteration, from the second iteration the ethr will be updated automatically.
    - *Default*: 0.01

    [back to top](#input-file)
- start_wfc<a id="start-wfc"></a>
    - *Type*: String
    - *Description*: Only useful for plane wave basis only now. It is the name of the starting wave functions. In the future we should also make this         variable available for localized orbitals set.
            - atomic:
            - file:
    - *Default*:atomic

    [back to top](#input-file)
- start_charge<a id="start-charge"></a>
    - *Type*: String
    - *Description*: This variable is used for both plane wave set and localized orbitals set. It indicates the type of starting density. If set this to ‘atomic’, the density is starting from summation of atomic density of single atoms. If set this to ‘file’, the density will be read in from file. The file should be in the output directory. Besides, when you do ‘nspin=1’ calculation, you only need the density file SPIN1_CHGCAR. However, if you do ‘nspin=2’ calculation, you also need the density file SPIN2_CHGCAR. The density file should be output with these names if you set out_charge = 1 in INPUT file.
        -  atomic:
        - file:
    - *Default*:atomic

    [back to top](#input-file)
### Electronic structure and geometry relaxation
This part of variables are used to control the electronic structure and geometry relaxation
calculations.
- basis_type<a id="basis-type"></a>
    - *Type*: String
    - *Description*: This is very important parameters to choose basis set in ABACUS.
        - *pw*: Using plane-wave basis set only.
        - *lcao_in_pw*: Expand the localized atomic set in plane-wave basis.
        - lcao: Using localized atomic orbital sets.
    - *Default*: pw

    [back to top](#input-file)
- ks_solver<a id="ks-solver"></a>
    - *Type*: String
    - *Description*: It’s about choice of diagonalization methods for hamiltonian matrix expanded in a certain basis set.

        For plane-wave basis,
        - cg: cg method.
        - david: david is the Davidson algorithm.

        For atomic orbitals basis,
        - genelpa: This method should be used if you choose localized orbitals.
        - hpseps: old method, still used.
        - lapack: lapack can be used for localized orbitals, but is only used for single processor.

        If you set ks_solver=‘hpseps’ for basis_type=‘pw’, the program will be stopped with an error message:

        ```
        hpseps can not be used with plane wave basis.
        ```

        Then the user has to correct the input file and restart the calculation.
    - *Default*: cg (pw) or genelpa (lcao)

    [back to top](#input-file)
- smearing<a id="smearing"></a>
    - *Type*: String
    - *Description*: It indicates which occupation and smearing method is used in the calculation.
        - fixed: use fixed occupations.
        - gauss or gaussian: use gaussian smearing method.
        - mp: use methfessel-paxton smearing method.
    - *Default*: fixed.

    [back to top](#input-file)
- sigma<a id="sigma"></a>
    - *Type*: Real
    - *Description*: energy range for smearing, the unit is Rydberg.
    - *Default*: 0.001

    [back to top](#input-file)
- mixing_type<a id="mixing-type"></a>
    - *Type*: String
    - *Description*: Charge mixing methods.
        - plain: Just simple mixing.
        - kerker: Use kerker method, which is the mixing method in G space.
        - pulay: Standard Pulay method.
        - pulay-kerker:
    - *Default*: pulay.

    [back to top](#input-file)
- mixing_beta<a id="mixing-beta"></a>
    - *Type*: Real
    - *Description*: mixing parameter: 0 means no new charge
    - *Default*: 0.7

    [back to top](#input-file)
- mixing_ndim<a id="mixing-ndim"></a>
    - *Type*: Integer
    - *Description*: It indicates the mixing dimensions in Pulay, Pulay method use the density from previous mixing_ndim steps and do a charge mixing based on these density.
    - *Default*: 8

    [back to top](#input-file)
- mixing_gg0<a id="mixing-gg0"></a>
    - *Type*: Real
    - *Description*: used in pulay-kerker mixing method
    - *Default*: 1.5

    [back to top](#input-file)
- gamma_only<a id="gamma-only"></a>
    - *Type*: Integer
    - *Description*: It is an important parameter in localized orbitals set.
    It you set gamma_only = 1, ABACUS use gamma only, the algorithm is fast and you don’t need to specify the k-points file. If you set gamma_only = 0, more than one k-point is used and the ABACUS is slower compared to gamma only algorithm.
    - *Default*: 0

    [back to top](#input-file)
- printe<a id="printe"></a>
    - *Type*: Integer
    - *Description*: Print out energy for each band for every printe steps
    - *Default*: 100

    [back to top](#input-file)
- niter<a id="niter"></a>
    - *Type*: Integer
    - *Description*:This variable indicates the maximal iteration number for electronic iterations.
    - *Default*: 40

    [back to top](#input-file)
- diago_cg_maxiter<a id="diago-cg-maxiter"></a>
    - *Type*: Integer
    - *Description*: Only useful when you use ks_solver = cg or ks_solver = david. It indicates the maximal iteration number for cg/david method.
    - *Default*: 40

    [back to top](#input-file)
- diago_david_ndim<a id="diago-david-ndim"></a>
    - *Type*: Integer
    - *Description*: Only useful when you use ks_solver = david. It indicates the maximal dimension for david method.
    - *Default*: 10

    [back to top](#input-file)
- dr2<a id="dr2"></a>
    - *Type*: Real
    - *Description*: An important parameter in ABACUS. It’s the threshold for electronic iteration. It represents the charge density error between two sequential density from electronic iterations. Usually for local orbitals, usually 1e-6 may be accurate enough.
    - *Default*:1e-06

    [back to top](#input-file)
- charge_extrap<a id="charge-extrap"></a>
    - *Type*: String
    - *Description*: Methods to do extrapolation of density when ABACUS is doing geometry relaxations.
        - atomic: atomic extrapolation
        - first-order: first-order extrapolation
        - second-order: second-order extrapolation
    - *Default*:atomic

    [back to top](#input-file)
- out_charge<a id="out-charge"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the charge density on real space grid. The name of the density file is SPIN1_CHGCAR and SPIN2_CHGCAR (if nspin = 2). Suppose each density on grid has coordinate (x; y; z). The circle order of the density on real space grid is: z is the outer loop, then y and finally x (x is moving fastest).
    - *Default*:0

    [back to top](#input-file)
- out_potential<a id="out-potential"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the local potential on real space grid. The name of the file is SPIN1_POT and SPIN2_POT (if nspin = 2). If set to 2, ABACUS will output the electrostatic potential on real space grid. The name of the file is ElecStaticPot and ElecStaticP ot_AV E (along the z-axis).
    - *Default*:0

    [back to top](#input-file)
- out_dm<a id="out-dm"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the density matrix of localized orbitals, only useful for localized orbitals set. The name of the output file is SPIN1_DM and SPIN2_DM in the output directory.
    - *Default*: 0

    [back to top](#input-file)
- out_wf<a id="out-wf"></a>
    - *Type*: Integer
    - *Description*: Only used in localized orbitals set. When set this variable to 1, it outputs the coefficients of wave functions. The file name is LOWF_GAMMA_S1:dat when using only GAMMA k-point for spin 1 component and LOWF_GAMMA_S2:dat for spin 2 component.
    - *Default*: 0

    [back to top](#input-file)
- out_lowf<a id="out-lowf"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the wave functions coefficients.
    - *Default*: 0

    [back to top](#input-file)
- mulliken<a id="mulliken"></a>
    - *Type*: Integer
    - *Description*: If set to 1, ABACUS will output the Mulliken population analysis result. The name of the output file is mulliken.txt
    - *Default*: 0

    [back to top](#input-file)
- out_alllog<a id="out-alllog"></a>
    - *Type*: Integer
    - *Description*: When set to 1, ABACUS will generate a log file for each processor when parallel, it is very useful for debugging.
    - *Default*: 0

    [back to top](#input-file)
- force<a id="force"></a>
    - *Description*: If set to 1, calculate the force at the end of the electronic iteration. 0 means the force calculation is turned off.
    - *Default*: 0

    [back to top](#input-file)
- nstep<a id="nstep"></a>
    - *Type*: Integer
    - *Description*: The maximal number of ionic iteration steps, the minimal value is 1.
    - *Default*: 1

    [back to top](#input-file)
- force_thr<a id="force-thr"></a>
    - *Type*: Real
    - *Description*: The threshold of the force convergence, it indicates the largest force among all the atoms, the unit is Ry=Bohr, *Default*: 0.000388935 Ry/Bohr = 0.01 eV/Angstrom

    [back to top](#input-file)
- force_thr_ev<a id="force-thr-ev"></a>
    - *Type*: Real
    - *Description*: The threshold of the force convergence, has the same function as force_thr, just the unit is different, it is eV=Angstrom, you can choose either one as you like. The recommendation value for using atomic orbitals is 0:04 eV/Angstrom.
    - *Default*: 0.01 eV/Angstrom

    [back to top](#input-file)
- bfgs_w1<a id="bfgs-w1"></a>
    - *Type*: Real
    - *Description*: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
    - *Default*: 0.01

    [back to top](#input-file)
- bfgs_w2<a id="bfgs-w2"></a>
    - *Type*: Real
    - *Description*: This variable controls the Wolfe condition for BFGS algorithm used in geometry relaxation. You can look into paper Phys.Chem.Chem.Phys.,2000,2,2177 for more information.
    - *Default*: 0.5

    [back to top](#input-file)
- trust_radius_max<a id="trust-radius-max"></a>
    - *Type*: Real
    - *Description*: This variable is for geometry optimization. It indicates the maximal movement of all the atoms. The sum of the movements from all atoms can be increased during the optimization steps. However, it will not be larger than trust_radius_max Bohr.
    - *Default*: 0.8

    [back to top](#input-file)
- trust_radius_min<a id="trust-radius-min"></a>
    - *Type*: Real
    - *Description*: This variable is for geometry optimization. It indicates the minimal movement of all the atoms. When the movement of all the atoms is smaller than trust_radius_min Bohr , and the force convergence is still not achieved, the calculation will break down.
    - *Default*: 1e-5

    [back to top](#input-file)
- trust_radius_ini<a id="trust-radius-ini"></a>
    - *Type*: Real
    - *Description*: This variable is for geometry optimization. It indicates the initial movement of all the atoms. The sum of the movements from all atoms is trust_radius_ini Bohr.
    - *Default*: 0.5

    [back to top](#input-file)
- stress<a id="stress"></a>
    - *Type*: Integer
    - *Description*: If set to 1, calculate the stress at the end of the electronic iteration. 0 means the stress calculation is turned off.
    - *Default*: 0

    [back to top](#input-file)
- stress_thr<a id="stress-thr"></a>
    - *Type*: Real
    - *Description*: The threshold of the stress convergence, it indicates the largest stress among all the directions, the unit is KBar,
    - *Default*: 10

    [back to top](#input-file)
- press1, 2, 3<a id="press"></a>
    - *Type*: Real
    - *Description*: the external pressures along three axes,the compressive stress is taken to be positive, the unit is KBar.
    - *Default*: 0

    [back to top](#input-file)
- fixed_axes<a id="fixed-axes"></a>
    - *Type*: String
    - *Description*:which axes are fixed when do cell relaxation. Possible choices are:
        - None : default; all can relax
        - volume : relaxation with fixed volume
        - a : fix a axis during relaxation
        - b : fix b axis during relaxation
        - c : fix c axis during relaxation
        - ab : fix both a and b axes during relaxation
        - ac : fix both a and c axes during relaxation
        - bc : fix both b and c axes during relaxation
        - abc : fix all three axes during relaxation
    - *Default*: None

    [back to top](#input-file)
- move_method<a id="move-method"></a>
    - *Type*: String
    - *Description*: The method to do geometry optimizations. If set to bfgs, using BFGS algorithm. If set to cg, using cg algorithm. If set to sd, using steepest-descent lgorithm.
    - *Default*: cg

    [back to top](#input-file)

### Molecular dynamics
This part of variables are used to control the molecular dynamics calculations.

- md_type<a id="md-type"></a>
    - *Type*: Integer
    - *Description*: control the ensemble to run md.
        - 0: When set to 0, ABACUS will use NVE ensemble;
        - 1: When set to 1, ABACUS will use NVT ensemble with Nose Hoover method;
        - 2: When set to 2, ABACUS will use NVT ensemble with Velosity Scaling method;
    - *Default*: 1

    [back to top](#input-file)
- md_rstmd<a id="md-rstmd"></a>
    - *Type*: Bool
    - *Description*: to control whether restart md.
        - 0:When set to 0, ABACUS will calculate md normolly.
        - 1:When set to 1, ABACUS will calculate md from last step in your test before.
    - *Default*: 0

    [back to top](#input-file)
- md_dt<a id="md_dt"></a>
    - *Type*: Double
    - *Description*: This is the time step(fs) used in md simulation .
    - *Default*: No default

    [back to top](#input-file)
- md_tfirst & md_tlast<a id="md-t"></a>
    - *Type*: Double
    - *Description*: This is the temperature used in md simulation, md_tlast’s default value is md_tfirst. If md_tlast is setted and be different from the md_tfirst, ABACUS will automatically generate a linear temperature gradient file named ”ChangeTemp.dat”, you can also set this file according to your needs instead.
    - *Default*: No default

    [back to top](#input-file)
- md_qmass<a id="md-qmass"></a>
    - *Type*: Double
    - *Description*: Inertia of extended system variable. Used only when md_type is 1 or 2, you should set a number which is larger than 0. If you want to autoset this by ABACUS,just set it to 0.
    - *Default*: 0

    [back to top](#input-file)
- md_nresn & md_nyosh<a id="md-nresn"></a>
    - *Type*: Integer
    - *Description*: Used when md_type is 1 or 2, control the Nose-Hoover thermostat extended-system, you can only set them at 1,3,5.
    - *Default*: md_nresn=md_nyosh=3

    [back to top](#input-file)
- md_dumpmdfred<a id="md_dumpmdfred"></a>
    - *Type*: Integer
    - *Description*:This is the steps to control the frequence to output md information
    - *Default*: 1

    [back to top](#input-file)
- md_domsd<a id="md_domsd"></a>
    - *Type*: Integer
    - *Description*:when set to 1, ABACUS will calculate mean square displacement and the diffusion of each element.
    - *Default*: 1

    [back to top](#input-file)
- md_fixtemperature<a id="md-fixtemperature"></a>
    - *Type*: Integer
    - *Description*:
        - n:when set to n (n > 1), ABACUS will read the file "ChangeTemp.dat" and change system’s temperature every n steps
        - 0,1:When set to 0 or 1, ABACUS won’t change the temperature during running MD.
    - *Default*: 1

    [back to top](#input-file)
- md_msdstarttime<a id="md-msdstarttime"></a>
    - *Type*: Integer
    - *Description*:when set to n, ABACUS will calculate mean square displacement and the diffusion of each element from nth step.
    - *Default*: 1

    [back to top](#input-file)

### VdW correction
This part of variables are used to control vdW-corrected related parameters.

- vdw_method<a id="vdw-method"></a>
    - *Type*: String
    - *Description*: If set to d2 ,d3_0 or d3_bj, ABACUS will calculate corresponding vdW correction, which is DFT-D2, DFT-D3(0) or DFTD3(BJ) method. And this correction includes energy and forces. "none" means that no vdW-corrected method has been used.
    - *Default*: none

    [back to top](#input-file)
- vdw_s6<a id="vdw-s6"></a>
    - *Type*: Real
    - *Description*: This scale factor is to optimize the interaction energy deviations. For DFT-D2, it is found to be 0.75 (PBE), 1.2 (BLYP), 1.05 (B-P86), 1.0 (TPSS), and 1.05 (B3LYP). For DFT-D3, recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: default

    [back to top](#input-file)
- vdw_s8<a id="vdw-s8"></a>
    - *Type*: Real
    - *Description*: This scale factor is only the parameter of DFTD3 approachs including D3(0) and D3(BJ). Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: default

    [back to top](#input-file)
-   vdw_a1<a id="vdw-a1"></a>
    - *Type*: Real
    - *Description*: This damping function parameter is only the parameter of DFT-D3 approachs including D3(0) and D3(BJ). Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: default

    [back to top](#input-file)
- vdw_a2<a id="vdw-a2"></a>
    - *Type*: Real
    - *Description*: This damping function arameter is only the parameter of DFT-D3(BJ) approach. Recommended values of this parameter with different DFT functionals can be found on the [webpage](https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3). The default values of this parameter in ABACUS is available for PBE. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: default

    [back to top](#input-file)
- vdw_d<a id="vdw-d"></a>
    - *Type*: Real
    - *Description*: The variable is to control the dumping speed of dumping function of DFT-D2.
    - *Default*: 20

    [back to top](#input-file)
- vdw_abc<a id="vdw-abc"></a>
    - *Type*: Integer
    - *Description*: The variable is to control the calculation of three-body term of DFT-D3 approachs, including D3(0) and D3(BJ). If set to 1, ABACUS will calculate three-body term, otherwise, the three-body term is not included.
    - *Default*: 0

    [back to top](#input-file)
- vdw_C6_file<a id="vdw-C6-file"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies the name of each elemetent’s C<sub>6</sub> Parameters file. If you don’t setup this, ABACUS will use the default C<sub>6</sub> Parameters stored in the programme already. We’ve stored elements from 1_H to 86_Rn. Otherwise, if you want to use some new C<sub>6</sub> Parameters, you should provide a file contained all the C<sub>6</sub> Parameters ordered by periodic table of elements, from 1_H to the last elements you want.
    - *Default*: default

    [back to top](#input-file)
- vdw_C6_unit<a id="vdw-C6-unit"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies unit of C<sub>6</sub> Parameters. Two kinds of unit is available: `J·nm6/mol` (means J·nm<sup>6</sup>/mol) and eVA(means eV·Angstrom)
    - *Default*: Jnm6/mol

    [back to top](#input-file)
- vdw_R0_file<a id="vdw-R0-file"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies the name of each elemetent’s R<sub>0</sub> Parameters file. If you don’t setup this, ABACUS will use the default R<sub>0</sub> Parameters stored in the programme already. We’ve stored elements from 1_H to 86_Rn. Otherwise, if you want to use some new R<sub>0</sub> Parameters, you should provide a file contained all the R<sub>0</sub> Parameters ordered by periodic table of elements, from 1_H to the last elements you want.
    - *Default*: default

    [back to top](#input-file)
- vdw_R0_unit<a id="vdw-R0-unit"></a>
    - *Type*: String
    - *Description*: This variable which is useful only when set vdw_method to d2 specifies unit of R<sub>0</sub> Parameters. Two kinds of unit is available: A(means Angstrom) and Bohr.
    - *Default*: A

    [back to top](#input-file)
- vdw_model<a id="vdw-model"></a>
    - *Type*: String
    - *Description*: To calculate the periodic structure, you can assign the number of lattice cells calculated. This variable specifies the kind of model to assign. If set to period, ABACUS will calculate a cubic periodic structure as assigned in vdw_period. If set to radius, ABACUS will calculate a cubic periodic structure containing a sphere, whose radius is vdw_radius and centre is origin point.
    - *Default*: radius

    [back to top](#input-file)
- vdw_radius<a id="vdw-radius"></a>
    - *Type*: Real
    - *Description*: If vdw_model is set to radius, this variable specifies the radius of the calculated sphere. For DFT-D2, the default value is 56.6918 ,while it is 95 for DFT-D3. This variable will be set to default, if no vdW-corrected method has been used.
    - *Default*: default

    [back to top](#input-file)
- vdw_radius_unit<a id="vdw-radius-unit"></a>
    - *Type*: String
    - *Description*: If vdw_model is set to radius, this variable specifies the unit of vdw_radius. Two kinds of unit is available: A(means Angstrom) and Bohr.
    - *Default*: Bohr

    [back to top](#input-file)
- vdw_cn_radius<a id="vdw-cn-radius"></a>
    - *Type*: Real
    - *Description*: This cutoff is chosen for the calculation of the coordination number (CN) in DFT-D3 approachs, including D3(0) and D3(BJ).
    - *Default*: 40.0

    [back to top](#input-file)
- vdw_cn_radius_unit<a id="vdw-cn-radius-unit"></a>
    - *Type*: String
    - *Description*: This variable specifies the unit of vdw_cn_radius. Two kinds of unit is available: A(means Angstrom) and Bohr.
    - *Default*: Bohr

    [back to top](#input-file)
- vdw_period<a id="vdw-period"></a>
    - *Type*: Int Int Int
    - *Description*: If vdw_model is set to period, these variables specify the number of x, y and z periodic.
    - *Default*: 3 3 3

    [back to top](#input-file)

### Berry phase and wannier90 interface
This part of variables are used to control berry phase and wannier90 interfacae parameters.

- berry_phase<a id="berry-phase"></a>
    - *Type*: Integer
    - *Description*: 1, calculate berry phase; 0, no calculate berry phase.
    - *Default*: 0

    [back to top](#input-file)
- gdir<a id="gdir"></a>
    - *Type*: Integer
    - *Description*:
        - 1: calculate the polarization in the direction of the lattice vector a_1 that is defined in STRU file.
        - 2: calculate the polarization in the direction of the lattice vector a_2 that is defined in STRU file.
        - 3: calculate the polarization in the direction of the lattice vector a_3 that is defined in STRU file.
    - *Default*: 3

    [back to top](#input-file)
- towannier90<a id="towannier90"></a>
    - *Type*: Integer
    - *Description*: 1, generate files for wannier90 code; 0, no generate.
    - *Default*: 0

    [back to top](#input-file)
- nnkpfile<a id="nnkpfile"></a>
    - *Type*: String
    - *Description*: the file name when you run “wannier90 -pp ...”.
    - *Default*: seedname.nnkp

    [back to top](#input-file)
- wannier_spin<a id="wannier-spin"></a>
    - *Type*: String
    - *Description*: If nspin is set to 2,
        - up: calculate spin up for wannier function.
        - down: calculate spin down for wannier function.
    - *Default*: up

    [back to top](#input-file)
- tddft<a id="tddft"></a>
    - *Type*: Integer
    - *Description*:
        - 1: calculate the real time time dependent density functional theory (TDDFT).
        - 0: do not calculate TDDFT.
    - *Default*: 0

    [back to top](#input-file)
- vext<a id="vext"></a>
    - *Type*: Integer
    - *Description*:
        - 1: add a laser material interaction (extern laser field).
        - 0: no extern laser field.
    - *Default*: 0

    [back to top](#input-file)
- vext_dire<a id="vext-dire"></a>
    - *Type*: Integer
    - *Description*:
        - 1: the direction of external light field is along x axis.
        - 2: the direction of external light field is along y axis.
        - 3: the direction of external light field is along z axis.
    - *Default*: 1

    [back to top](#input-file)