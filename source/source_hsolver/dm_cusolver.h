#ifndef ABACUS_DM_CUSOLVER_H
#define ABACUS_DM_CUSOLVER_H

#include "hs_matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"
#include <complex>
#include <string>

namespace hsolver
{
// Immutable inputs for the cuSOLVER density-matrix path.
struct CusolverDmConfig
{
    const std::string& method;
    const int kpar;
    const int nspin;
    const bool skip_charge;
    const int world_rank;
    const int bands;
};

// Return false for configurations using the existing CPU DM path; GPU failures are fatal.
bool solve_with_gpu_dm(HSMatrix<double>& hs, psi::Psi<double>& wave, elecstate::ElecState* state,
                       module_dm::DensityMatrix<double, double>& density,
                       const Parallel_Orbitals* orbitals, const CusolverDmConfig& config);
bool solve_with_gpu_dm(HSMatrix<std::complex<double>>& hs, psi::Psi<std::complex<double>>& wave,
                       elecstate::ElecState* state, module_dm::DensityMatrix<std::complex<double>, double>& density,
                       const Parallel_Orbitals* orbitals, const CusolverDmConfig& config);
}

#endif
