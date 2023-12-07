# Quasiatomic Orbital Analysis for ABACUS

- Hamiltonian matrix transformation

## Requirements

- [NumPy](https://numpy.org/)
- [SciPy](https://www.scipy.org/)
- [Matplotlib](https://matplotlib.org/)

## Installation

This code runs like scripts, so no need to install.

## Usage

1. Perform ABACUS `basis_type lcao` calculation with additional keywords:
   ```Plain text
   qo_switch 1
   qo_basis pswfc
   ```
   , this will tell ABACUS to extract pseudowavefunction from pseudopotential file. To use this without any errors, please make sure there IS pseudowavefunction in pseudopotential file.  
   One can also add `qo_screening_coeff` keyword in `INPUT` to tune the behavior of pseudowavefunction:
   ```Plain text
   qo_screening_coeff 0.1
   ```
   `qo_screening_coeff` always to be `0.1`, but whatever value you like, it must be a positive number to ensure a proper behavior at least at infinity, i.e. $\lim_{|\mathbf{r}|\rightarrow\infty}\tilde{\psi}(|\mathbf{r}|)\rightarrow0$.  
   There is also one *very experimental, very numerically instable* feature:
   ```Plain text
   qo_switch 1
   qo_basis hydrogen
   ```
   To use hydrogen-like radial function as projector. However, because the charges cannot be set very physically presently (unless let user set them), ABACUS will read atomic charge from pseudopotential file, which may cause unphysical shrink or expansion.
2. Copy output Hamiltonian and overlap matrices $H(\mathbf{k})$ and $S(\mathbf{k})$ files (`data-[i]-H` and `data-[i]-S`), the overlap matrices between AO and numerical atomic orbitals $S^{\chi\phi}(\mathbf{k})$ `QO_ovlp_[i].dat`, converged wavefunction `LOWF_K_[i].txt` from ABACUS output `OUT.*` directory to the path you like.
3. Specify the path, number of kpoints, number of states want to introduce and the lattice vector $\mathbf{R}$ (presently it is the one set to $(0, 0, 0)$) in `source/main.py`.
    ```python
    import components.driver as driver

    if __name__ == "__main__":

        path = "./examples/input_pswfc_si2"
        nkpts = 8
        band_range = (0, 13) # min and max indices of band you want to reproduce

        d_ = driver.toQO_Driver() # initiate the driver of toQO python-end
        d_.initialize(path, nkpts, band_range) # initialize with parameters set
        d_.space_expansion() # expand the space spanned by eigenstates
        d_.reproduce_hamiltonian() # reproduce the energy spectrum
    ```
4.  run it!
