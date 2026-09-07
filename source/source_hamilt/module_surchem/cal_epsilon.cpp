#include "surchem.h"

#include <cassert>

void surchem::cal_epsilon(const ModulePW::PW_Basis* rho_basis,
                          const double* PS_TOTN_real,
                          double* epsilon,
                          double* epsilon0)
{
    assert(this->parameters_set_
           && "surchem::set_parameters() must be called before using the solvent model");
    // shapefunction value varies from 0 in the solute to 1 in the solvent
    // epsilon = 1.0 + (eb_k - 1) * shape function
    // build epsilon in real space (nrxx)
    double *shapefunc = new double[rho_basis->nrxx];
    for (int i = 0; i < rho_basis->nrxx; i++)
    {
        shapefunc[i] = erfc((log(std::max(PS_TOTN_real[i], 1e-10) / this->parameters_.nc_k))
                            / sqrt(2.0) / this->parameters_.sigma_k) / 2;
        epsilon[i] = 1 + (this->parameters_.eb_k - 1) * shapefunc[i];
        epsilon0[i] = 1.0;
    }
    delete[] shapefunc;
}