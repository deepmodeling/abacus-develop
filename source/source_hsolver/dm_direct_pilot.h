#ifndef ABACUS_DM_DIRECT_PILOT_H
#define ABACUS_DM_DIRECT_PILOT_H

#include "hs_matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/elecstate.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"
#include <chrono>
#include <complex>
#include <cstddef>
#include <string>

namespace hsolver
{
// Immutable per-call selection, restricted to the private trial build.
struct DirectConfig
{
    const std::string& method;
    const int kpar;
    const int nspin;
    const bool skip_charge;
    const int world_rank;
    const int world_size;
    const int bands;
    const char* mode;
};

void direct_sync(int context);
double direct_elapsed(const std::chrono::steady_clock::time_point& started);
bool direct_gamma(HSMatrix<double>& hs, psi::Psi<double>& wave, elecstate::ElecState* state,
                  module_dm::DensityMatrix<double, double>& density,
                  const Parallel_Orbitals* orbitals, const DirectConfig& config);
bool direct_gamma(HSMatrix<std::complex<double>>& hs, psi::Psi<std::complex<double>>& wave,
                  elecstate::ElecState* state, module_dm::DensityMatrix<std::complex<double>, double>& density,
                  const Parallel_Orbitals* orbitals, const DirectConfig& config);
void direct_report_fallback(bool root, const char* mode);
void direct_report_cpu(bool root, int context, const std::chrono::steady_clock::time_point& started);
void direct_report_diagonal_start(bool root, const std::chrono::steady_clock::time_point& started);
void direct_report_redistribution(bool root, const std::chrono::steady_clock::time_point& started);
void direct_report_combined(bool root, int context, const char* mode, bool direct,
                            const std::chrono::steady_clock::time_point& started);
void direct_report_phase(int rank, int rows, int bands, std::size_t scalar_bytes,
                          double gather_seconds, double diagonal_seconds, double psi_seconds);
}

#endif
