"""Main program of Quasiatomic Orbital (QO) analysis postprocessing package

Author: ykhuang
Date: 2023/12/25, Merry Christmas

Instruction
-----------
This program is used to postprocess the output of ABACUS calculation. To
obtain files needed, please set the following parameters in ABACUS input file:
```
# Quasiatomic Orbital Analysis
qo_switch 1            # turn on quasiatomic orbital analysis
qo_basis pswfc         # use pseudowavefunction as QO corresponding AO
qo_thr 1e-6            # controls the realspace spreading of AO
qo_screening_coeff 0.5 # controls the exponential decay-manner of AO
```
"""
import abacus2qo.components.driver as driver
import numpy as np

if __name__ == "__main__":
    # example of QO
    path = "./examples/OUT.ABACUS/"
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