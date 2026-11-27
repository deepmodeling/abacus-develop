#ifndef ABACUS_DM_DIRECT_PILOT_H
#define ABACUS_DM_DIRECT_PILOT_H

#include "hs_matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"
#include <chrono>
#include <complex>
#include <string>

namespace hsolver
{
// Immutable selection passed explicitly from the INPUT-configured solve.
struct DirectConfig
{
    const std::string& method;
    const int kpar;
    const int nspin;
    const bool skip_charge;
    const int world_rank;
    const int bands;
};

void direct_sync(int context);
double direct_elapsed(const std::chrono::steady_clock::time_point& started);
void direct_gamma(HSMatrix<double>& hs, psi::Psi<double>& wave, elecstate::ElecState* state,
                  module_dm::DensityMatrix<double, double>& density,
                  const Parallel_Orbitals* orbitals, const DirectConfig& config);
void direct_gamma(HSMatrix<std::complex<double>>& hs, psi::Psi<std::complex<double>>& wave,
                  elecstate::ElecState* state, module_dm::DensityMatrix<std::complex<double>, double>& density,
                  const Parallel_Orbitals* orbitals, const DirectConfig& config);
void direct_report_diagonal_start(bool root, const std::chrono::steady_clock::time_point& started);
void direct_report_redistribution(bool root, const std::chrono::steady_clock::time_point& started);
}

#endif
