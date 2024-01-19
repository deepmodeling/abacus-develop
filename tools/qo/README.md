# Quasiatomic Orbital Analysis for ABACUS

- Hamiltonian matrix transformation

## Requirements

- [NumPy](https://numpy.org/)
- [SciPy](https://www.scipy.org/)
- [Matplotlib](https://matplotlib.org/) (optional, for plotting)

## Installation

This code runs like scripts, you can also install this by:
```bash
python3 setup.py install
```
, then you can use it like a package:
```python
import abacus2qo
```
If you would like to change code and run it, you can also install it by:
```bash
python3 setup.py develop
```

## Usage

1. Perform ABACUS `basis_type lcao` calculation with additional keywords:
   ```Plain text
   qo_switch 1
   qo_basis pswfc
   ```
   , this will tell ABACUS to extract pseudowavefunction from pseudopotential file. To use this without any errors, please make sure there IS pseudowavefunction in pseudopotential file.  
   One can also add `qo_screening_coeff` keyword in `INPUT` to tune the behavior of pseudowavefunction (default is set to 0.1):
   ```Plain text
   qo_screening_coeff 0.1
   ```
   `qo_screening_coeff` always to be `0.1`, but whatever value you like, it must be a positive number to ensure a proper behavior at least at infinity, i.e. $\lim_{|\mathbf{r}|\rightarrow\infty}\tilde{\psi}(|\mathbf{r}|)\rightarrow0$.  
   There is also a parameter `qo_thr` controlling the tailing of pseudowavefunction, which defines as the convergence threshold of norm of pseudowavefunction when realspace cutoff increasing. To change this parameter, just add it in `INPUT`:
   ```Plain text
   qo_thr 1e-6
   ```
   There is also one *very experimental, very numerically instable* feature:
   ```Plain text
   qo_switch 1
   qo_basis hydrogen
   ```
   To use hydrogen-like radial function as projector. However, because the charges cannot be set very physically presently (unless let user set them), ABACUS will read atomic charge from pseudopotential file, which may cause unphysical shrink or expansion.
   For example see `examples`, find `INPUT`, `STRU` and `KPT`.
2. Copy output Hamiltonian and overlap matrices $H(\mathbf{k})$ and $S(\mathbf{k})$ files (`data-[i]-H` and `data-[i]-S`), the overlap matrices between AO and numerical atomic orbitals $S^{\chi\phi}(\mathbf{k})$ `QO_ovlp_[i].dat`, converged wavefunction `LOWF_K_[i].txt`, kpoints information summary `kpoints` from ABACUS output `OUT.*` directory to the path you like. Or in another way, just assign the `path` in `source/main.py` to the path of `OUT.*` directory.
3. Then you specify the path, number of kpoints, number of states want to introduce and the lattice vector $\mathbf{R}$ (presently it is the one set to $(0, 0, 0)$) in `source/main.py`.
    ```python
    import source.components.driver as driver
    import numpy as np

    if __name__ == "__main__":
        # example of QO
        path = "./tests/integrate/220_NO_KP_QO/OUT.ABACUS"
        # fully non-reduced kpoints
        nkpts = 125
        # band range to reproduce
        band_range = (0, 13)
        
        # create driver
        d_ = driver.toQO_Driver()
        # initialize driver
        d_.initialize(path, nkpts, "scf", band_range)
        # expand space spanned by selected bands
        d_.space_expansion()
        # calculate H(R) and S(R) in R-space in QO representation
        HRs, SRs = d_.reproduce_hamiltonian(Rs=[
            np.array([0, 0, 0])
            ])
    ```
    `HRs` and `SRs` are list of `numpy.ndarray`, so it is convinient to use `numpy.save` to save them.
    You can also find a file `QO_supercells.dat` in `OUT.${suffix}/`, which contains the information of supercell. You can use it as a reference to construct `Rs`.
4.  run it!

## Formulation info

QO is simple in formulation, its central idea is to maximize the norm of QO, which is the projection of AO in extended space. The space is spanned by the selected bands concatenated with the pseudowavefunction.  

First project pseudowavefunction with NAO, which means add projection operator of NAO at left side of AO:  
$\hat{P}|A_\mu(\mathbf{k})\rangle$  
, where $\mathbf{k}$ specifies the kpoint, $\mu$ distinguishes the AO, it is actually a index describing which type, which atom, which $(n, l, m)$. Projection operator $\hat{P}$ can be represented as:  
$\hat{P} = \sum_{\alpha\beta}{|\phi_{\alpha}(\mathbf{k})\rangle \mathbf{S}^{-1}_{\alpha\beta}\langle\phi_{\beta}(\mathbf{k})|A_{\mu}(\mathbf{k})\rangle}$  
The matrix element $\mathbf{S}^{-1}_{\alpha\beta}$, hopefully denoted as it is here without any misleading, is the matrix element of inverse overlap matrix between NAOs.  

Similarly, we project AO into selected eigenstates spanned space:  
$\hat{P}_{\psi}|A_\mu(\mathbf{k})\rangle$ = $\sum_{n}{|\psi_{\alpha}(\mathbf{k})\rangle\langle\psi_{\beta}(\mathbf{k})|}A_{\mu}(\mathbf{k})\rangle$  
, one can quickly rewrite it as:  
$\begin{cases}
	|A_{\mu}(\mathbf{k})\rangle \rightarrow \sum_{\alpha \beta}{\left( \mathbf{S}_{\alpha \beta}^{-1}\langle \phi _{\beta}(\mathbf{k})|A_{\mu}(\mathbf{k})\rangle \right) |\phi _{\alpha}(\mathbf{k})\rangle}\\
	|A_{\mu}^{\parallel}(\mathbf{k})\rangle \rightarrow \sum_{\alpha}{\left( \sum_{\beta}{\langle \phi _{\beta}(\mathbf{k})|A_{\mu}(\mathbf{k})\rangle \sum_n{c_{n\alpha}\left( \mathbf{k} \right) c_{n\beta}^{\dagger}\left( \mathbf{k} \right)}} \right) |\phi _{\alpha}(\mathbf{k})\rangle}\\
\end{cases}$  

Therefore the virtual states introduced by AO, which definitely is NAOs linearly combined in some manner, is:  
$|A_{\mu}^{\bot}(\mathbf{k})\rangle =|A_{\mu}(\mathbf{k})\rangle -|A_{\mu}^{\parallel}(\mathbf{k})\rangle$  
Then canonical orhogonalization is performed, get only parts of eigenstates, denoted as $\{c_m(\mathbf{k})\}$. Then merge spaces spanned by $\{c_m(\mathbf{k})\}$ and $\{|\psi_{n}(\mathbf{k})\rangle\}$, denote as $\{\varphi_i(\mathbf{k})\}$.  
Then the QO is:  
$|Q_{\mu}(\mathbf{k})\rangle =\sum_i{|\varphi _i(\mathbf{k})\rangle \langle \varphi _i(\mathbf{k})|A_{\mu}(\mathbf{k})\rangle}$