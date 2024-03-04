*Use Markdown render to convert this document into directly human-readable version*
# ABACUS LCAO2QO Hands-on
## Introduction
This document is a hands-on tutorial for the ABACUS LCAO2QO tool. For basic information about QO (Quasiatomic Orbital) method, please refer to *Quasiatomic orbitals for ab initio tight-binding analysis* by *Qian et al.*: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.245112.
### Basic ideas
QO is for projecting as many as possible eigenstates of a system onto a atom-centered basis set, unlike MLWFs (Maximally Localized Wannier Functions) that would produce bond-like Wannier functions for some systems. The QO is defined from an optimization problem that maximal the norm of the projection of atom-centered orbitals onto the super-plane spanned by eigenstates and constructed unoccupied states:  
$\max\mathcal{L}:=\max \sum_{\mu}{\left\|\left(\sum_{n\mathbf{k}}{\hat{P}_{\psi_{n\mathbf{k}}}}+\sum_{m\mathbf{k}}{\hat{P}_{c_{m\mathbf{k}}}}\right)|A_{\mu}\rangle\right\|^2}$  
where $\hat{P}_{\psi_{n\mathbf{k}}}$ and $\hat{P}_{c_{m\mathbf{k}}}$ are the projection operators for occupied and unoccupied states, respectively. $A_{\mu}$ is the atom-centered orbital.  
The former part of projection is fixed once the occupied states are determined, and the latter part is determined by selection of $Nq-Nb$ from $\sum_{m\mathbf{k}}{\hat{P}_{c_{m\mathbf{k}}}}|A_{\mu}\rangle$ states, with total number of $Nq$ and number of occupied states $Nb$.  
### Formulation
See Github issue #3640 for details: https://github.com/deepmodeling/abacus-develop/issues/3640  
## Prerequisites
- ABACUS
- Python3
- Numpy
- Scipy
- Matplotlib (optional)
## Usage
### Basis set selection
ABACUS provides various kinds of atom-centered orbitals construction, such as hydrogen-like orbitals, exponentially decay scaled pseudowavefunctions, single zeta (not implemented, develop up to user actual needs). Select orbitals that fit most to the system you are interested in.
### Input file of ABACUS
To run QO analysis, in ABACUS INPUT file, following keywords are needed to be set:
```bash
qo_switch           1
qo_basis            hydrogen
qo_strategy         energy-valence
qo_thr              1.0e-10
#qo_screening_coeff 0.1
```
- `qo_switch` is a switch to turn on QO analysis.
- `qo_basis` is the basis set selection, parameters above define the hydrogen-like orbital construction
- `qo_strategy` is the strategy to set what orbitals needed to construct. For hydrogen-like orbitals, for example if element is Cu (1s2 2s2 2p6 3s2 3p6 3d10 4s1), `energy-valence` corresponds to construct only 3p 3d and 4s orbitals, `energy-full` corresponds to construct all orbitals, `minimal-nodeless` corresponds to construct 1s, 2p, 3d and 4f orbitals, `minimal-valence` corresponds to construct 4s 4p 4d 4f orbitals. For more information about this keyword, please refer to ABACUS manual: https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html#qo-strategy
- `qo_thr` is the threshold to determine the realspace tailing of orbitals, which is used to determine the cutoff radius of orbitals. Default value is 1.0e-10.
- `qo_screening_coeff`, for `qo_basis hydrogen` it controls whether use Slater screening coefficients to include many-electron effect into the orbitals. Any float number as value will turn on this feature. Default value is off.  

For some systems, pseudowavefunction works better than hydrogen-like orbitals (but not all pseudopotential has pseudowavefunction in file, check your pseudopotential file before run). To use pseudowavefunction as `qo_basis`, set as following:
```bash
qo_switch           1
qo_basis            pswfc
qo_strategy         all
qo_thr              1.0e-10
qo_screening_coeff  0.5
```
For `qo_basis pswfc`, `qo_screening_coeff` is necessary and for Si single crystal the optimal value is tested to be `0.5`. Test on this parameter should be done for different systems. A larger value will result in more localized orbitals.
### Output file of ABACUS
Once the calculation is done, the output file will contain the QO analysis results. The QO analysis results are stored in `OUT.${suffix}` directory, where `${suffix}` is set in `INPUT` and default is `ABACUS`. You will find following files:
```bash
data-0-H ...
data-0-S ...
```
These files are Hamiltonian and overlap matrix of numerical atomic orbitals basis sets. For each kpoint a pair of these files will be generated, and different spin channels are defined as different kpoints.    
```bash
LOWF_K_1.txt ...
```
These files are eigenstates represented by numerical atomic orbital basis sets. For each kpoint a file will be generated.  
```bash
QO_ovlp_0.dat ...
```
This file contains information of $\langle A_\mu(\mathbf{k})|\phi_\nu(\mathbf{k})\rangle$, file name indiced by kpoint.  
```bash
QO_supercells.dat
```
This file stores information of supercell considered when calculating two-center-integrals, will be usefull in kpoint extrapolation.  
```bash
kpoints
istate.info
running_scf.log
...
```
Basic output files of ABACUS.
### Run `postprocess.py`
After the calculation is done, run `postprocess.py` to generate QO analysis results. The script will read the output files of ABACUS and generate the QO analysis results. Following defines an example of `postprocess.py`:
```python
if __name__ == "__main__":

    # False to run unittest, True to run production
    run_type = "production" # can be "unittest" or "production"
    ndim = 7 # dimension of Monkhorst-Pack grid
    qo_basis = "hydrogen"
    qo_strategy = "minimal-valence"
    path = f"./scf/{qo_basis}/{qo_strategy}/OUT.ABACUS-{ndim}{ndim}{ndim}/"
    #fpic = f"{qo_basis}-{qo_strategy}-{ndim}{ndim}{ndim}.png"
    # only for production
    eig_thr = 1e-10
    ib_min, ib_max = 0, 4

    if run_type == "production":
        nkpts = ndim*ndim*ndim
        hqos_k, sqos_k, kvecs_d = production(path, nkpts, ib_min, ib_max, eig_thr, error_estimation_with_lcao=True)

    elif run_type == "unittest":
        unittest.main()
```
The `production` function will read the output files of ABACUS and generate the QO analysis results. The `path` is the path to the output files of ABACUS. `nkpts` is the number of kpoints. `ib_min` and `ib_max` are the range of orbitals to be considered. `eig_thr` is the threshold to determine the eigenstates to be considered. `error_estimation_with_lcao` is a switch to turn on error estimation with LCAO. If set to `True`, the script will calculate the error estimation with LCAO. Collect information output by `production(...)` function. The `hqos_k` and `sqos_k` are the Hamiltonian and overlap matrix of QO basis sets. The `kvecs_d` is the direct coordinates of kpoints.  
QO postprocess is programmed with unittest, therefore to run the script, set `run_type` to `production` and run the script (but mostly you don't need this).
## Acknowledgement
ABACUS LCAO2QO module uses two-center-integrator module refactored by @jinzx10, who also helps in code debugging. @mohanchen, @dyzheng, @WHUweiqingzhou and @QG-phys provide suggestions on code maintainance and participate discussion in technical aspects.