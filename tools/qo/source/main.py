import components.driver as driver

if __name__ == "__main__":

    path = "./examples/input_pswfc_si2"
    nkpts = 8
    band_range = (0, 13)

    d_ = driver.toQO_Driver()
    d_.initialize(path, nkpts, band_range)
    d_.space_expansion()
    d_.reproduce_hamiltonian()