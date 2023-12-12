import source.components.driver as driver
import numpy as np

if __name__ == "__main__":

    path = "./examples/pswfc_scf_lcao_ZnO"
    nkpts = 30
    band_range = (0, 13)

    d_ = driver.toQO_Driver()
    d_.initialize(path, nkpts, "scf", band_range)
    d_.space_expansion()
    d_.reproduce_hamiltonian(Rs=[
        np.array([0, 0, 0])
        ])