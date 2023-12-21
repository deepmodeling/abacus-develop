import source.components.driver as driver
import numpy as np

if __name__ == "__main__":

    path = "/root/abacus-develop/dev/tests/integrate/220_NO_KP_QO/OUT.ABACUS"
    nkpts = 125
    band_range = (0, 13)

    d_ = driver.toQO_Driver()
    d_.initialize(path, nkpts, "scf", band_range)
    d_.space_expansion()
    d_.reproduce_hamiltonian(Rs=[
        np.array([0, 0, 0])
        ])