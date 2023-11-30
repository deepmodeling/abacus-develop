# Quasiatomic Orbital Analysis for ABACUS

- Hamiltonian matrix transformation

## Requirements

- [NumPy](https://numpy.org/)
- [SciPy](https://www.scipy.org/)
- [Matplotlib](https://matplotlib.org/)

## Installation

This code runs like scripts, so no need to install.

## Usage

1. Perform ABACUS `basis_type lcao` calculation with keyword `qo_switch 1` (optionally `qo_basis hydrogen` is for setting hydrogen-like basis functions, `qo_thr` is for searching the realspace cutoff radius of those orbitals, default value is `1e-6`).
2. Copy output Hamiltonian and overlap matrices $H(\mathbf{k})$ and $S(\mathbf{k})$ files (`data-[i]-H` and `data-[i]-S`), the overlap matrices between AO and numerical atomic orbitals $S^{\chi\phi}(\mathbf{k})$ `QO_ovlp_[i].dat`, converged wavefunction `LOWF_K_[i].txt` from ABACUS output `OUT.*` directory to the path you like.
3. Specify the path, number of kpoints, number of states want to introduce and the lattice vector $\mathbf{R}$ (presently it is the one set to $(0, 0, 0)$) in `source/to_qo.py`.
    ```python
    if __name__ == "__main__":
        # user settings
        nkpts = 8 # number of kpoints in ABACUS calculation
        m = 2 # number of new states added
        supercell = np.array([0, 0, 0]) # supercell coordinate

        tqo = toQO()
        tqo.initialize(nkpts, m, "./examples/input/")
        # calculate
    ```
4.  and run it! You will get output files like `hqoR0_0_0.txt`. 
