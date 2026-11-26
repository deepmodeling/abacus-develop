#include "dm_direct_pilot.h"
#include "kernels/cuda/diag_cusolver.cuh"
#include "source_base/module_device/device.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/tool_quit.h"
#include "source_estate/elecstate_tools.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

namespace hsolver
{
namespace
{
bool direct_supported(const int context, const int kpoints, const DirectConfig& config)
{
    const base_device::DeviceContext& device_context = base_device::DeviceContext::instance();
    int unsupported_rank = device_context.get_device_count() != 1
                   || device_context.get_local_rank() != config.world_rank;
    char scope = 'A';
    char topology = ' ';
    Cigsum2d(context, &scope, &topology, 1, 1, &unsupported_rank, 1, -1, -1);
    return config.method == "cusolver" && config.kpar == 1
                   && config.nspin == 1 && kpoints == 1 && !config.skip_charge
                   && unsupported_rank == 0;
}
}

void direct_sync(const int context)
{
    int token = 1;
    char scope = 'A';
    char topology = ' ';
    Cigsum2d(context, &scope, &topology, 1, 1, &token, 1, -1, -1);
}

double direct_elapsed(const std::chrono::steady_clock::time_point& started)
{
    const std::chrono::steady_clock::time_point finished = std::chrono::steady_clock::now();
    const std::chrono::steady_clock::duration duration = finished - started;
    return std::chrono::duration<double>(duration).count();
}

void direct_gamma(HSMatrix<double>& hs, psi::Psi<double>& wave, elecstate::ElecState* state,
                  module_dm::DensityMatrix<double, double>& density, const Parallel_Orbitals* orbitals,
                  const DirectConfig& config)
{
    const int bands = config.bands;
    const Parallel_Blkcyclic& layout = orbitals->para_2d();
    const int context = layout.ictxt();
    const int kpoints = wave.get_nk();
    const bool supported = direct_supported(context, kpoints, config);
    if (!supported)
    {
        ModuleBase::WARNING_QUIT("direct_gamma",
                                "device_dm gpu requires real Gamma FP64 LCAO with ks_solver cusolver, "
                                "nspin 1, kpar 1, charge updates, MPI ranks on one node and "
                                "one visible GPU per rank; select device_dm cpu in INPUT otherwise");
    }
    const int n = layout.nrglb();
    const bool root = layout.irproc() == 0 && layout.icproc() == 0;
    const std::size_t matrix_size = static_cast<std::size_t>(n) * n;
    const std::size_t wave_size = static_cast<std::size_t>(n) * bands;
    const std::size_t matrix_bytes = matrix_size * sizeof(double);
    const std::size_t hs_bytes = 2 * matrix_bytes;
    const std::size_t wave_bytes = wave_size * sizeof(double);
    const std::size_t root_matrix_size = root ? matrix_size : 1;
    const std::size_t root_wave_size = root ? wave_size : 1;
    std::vector<double> h(root_matrix_size);
    std::vector<double> s(root_matrix_size);
    std::vector<double> host_wave;
    std::vector<double> eigenvalues(n);
    std::vector<double> host_dm;
    const int leading_dimension = layout.irproc() == 0 ? n : 1;
    int full_desc[9] = {1, context, n, n, n, n, 0, 0, leading_dimension};
    int wave_desc[9] = {1, context, n, bands, n, bands, 0, 0, leading_dimension};
    int distributed_desc[9];
    int distributed_wave_desc[9];
    const int* matrix_descriptor = layout.desc();
    const int* matrix_descriptor_end = matrix_descriptor + 9;
    const int* wave_descriptor = orbitals->desc_wfc;
    const int* wave_descriptor_end = wave_descriptor + 9;
    std::copy(matrix_descriptor, matrix_descriptor_end, distributed_desc);
    std::copy(wave_descriptor, wave_descriptor_end, distributed_wave_desc);
    const std::chrono::steady_clock::time_point gather_start = std::chrono::steady_clock::now();
    ModuleBase::MatrixBlock<double> hk;
    ModuleBase::MatrixBlock<double> sk;
    hs.hs_at_k(0, hk, sk);
    direct_sync(context);
    const std::chrono::steady_clock::time_point redistribution_start = std::chrono::steady_clock::now();
    Cpdgemr2d(n, n, hk.p, 1, 1, distributed_desc, h.data(), 1, 1, full_desc, context);
    Cpdgemr2d(n, n, sk.p, 1, 1, distributed_desc, s.data(), 1, 1, full_desc, context);
    direct_sync(context);
    direct_report_redistribution(root, redistribution_start);
    const double gather_seconds = direct_elapsed(gather_start);
    const std::chrono::steady_clock::time_point diagonal_start = std::chrono::steady_clock::now();
    direct_report_diagonal_start(root, diagonal_start);
    host_wave.resize(root_wave_size);
    std::unique_ptr<Diag_Cusolver_gvd> solver;
    if (root)
    {
        solver.reset(new Diag_Cusolver_gvd());
        bool solved = false;
        const bool subset = Diag_Cusolver_gvd::use_subset(n, bands);
        if (subset)
        {
            solved = solver->solve_subset(n, bands, h.data(), s.data(), eigenvalues.data(), host_wave.data(), true);
        }
        if (!solved)
        {
            solver->solve_full(n, bands, h.data(), s.data(), eigenvalues.data(), host_wave.data(), true);
        }
    }
    direct_sync(context);
    const double diagonal_seconds = direct_elapsed(diagonal_start);
    const std::chrono::steady_clock::time_point psi_start = std::chrono::steady_clock::now();
    char scope = 'A';
    char topology = ' ';
    if (root) Cdgebs2d(context, &scope, &topology, bands, 1, eigenvalues.data(), bands);
    else Cdgebr2d(context, &scope, &topology, bands, 1, eigenvalues.data(), bands, 0, 0);
    for (int band = 0; band < bands; ++band) state->ekb(0, band) = eigenvalues[band];
    wave.fix_k(0);
    double* wave_pointer = wave.get_pointer();
    Cpdgemr2d(n, bands, host_wave.data(), 1, 1, wave_desc,
              wave_pointer, 1, 1, distributed_wave_desc, context);
    direct_sync(context);
    const double psi_seconds = direct_elapsed(psi_start);
    // The owning solver remains alive through ABACUS's actual occupation code.
    const std::chrono::steady_clock::time_point occupation_start = std::chrono::steady_clock::now();
    elecstate::calculate_weights(state->ekb, state->wg, state->klist, state->eferm,
                                 state->f_en, state->nelec_spin, bands, state->skip_weights);
    elecstate::calEBand(state->ekb, state->wg, state->f_en);
    direct_sync(context);
    const double occupation_seconds = direct_elapsed(occupation_start);
    const std::chrono::steady_clock::time_point dm_start = std::chrono::steady_clock::now();
    host_dm.resize(root_matrix_size);
    int success = 1;
    std::string error;
    if (root)
    {
        const double* weights = &state->wg(0, 0);
        success = solver->density_from_retained(bands, weights, host_dm.data(), error);
    }
    if (root) Cigebs2d(context, &scope, &topology, 1, 1, &success, 1);
    else Cigebr2d(context, &scope, &topology, 1, 1, &success, 1, 0, 0);
    if (!success)
    {
        const MPI_Comm comm = layout.comm();
        const Parallel::ParaWorld world = Parallel::ParaWorld::make_mpi("direct_dm", comm);
        Parallel::bcast_string(error, world, 0);
        ModuleBase::WARNING_QUIT("direct_gamma", error, EXIT_FAILURE);
    }
    const double device_dm_seconds = direct_elapsed(dm_start);
    const std::chrono::steady_clock::time_point scatter_start = std::chrono::steady_clock::now();
    double* dm_pointer = density.get_dmk_ptr(0);
    Cpdgemr2d(n, n, host_dm.data(), 1, 1, full_desc,
              dm_pointer, 1, 1, distributed_desc, context);
    direct_sync(context);
    const double scatter_seconds = direct_elapsed(scatter_start);
    if (root)
    {
    std::printf("DM_DIRECT_PHASE mode=direct n=%d bands=%d gather_hs=%.9e diagonal=%.9e psi_redistribute=%.9e "
                    "occupations=%.9e device_dm=%.9e dm_redistribute=%.9e hs_redistributions=2 hs_logical_bytes=%zu "
                    "psi_redistributions=1 psi_logical_bytes=%zu dm_redistributions=%d dm_logical_bytes=%zu\n",
                    n, bands, gather_seconds, diagonal_seconds, psi_seconds, occupation_seconds,
                    device_dm_seconds, scatter_seconds, hs_bytes, wave_bytes, success, matrix_bytes);
    }
    return;
}

void direct_gamma(HSMatrix<std::complex<double>>&, psi::Psi<std::complex<double>>&,
                  elecstate::ElecState*, module_dm::DensityMatrix<std::complex<double>, double>&,
                  const Parallel_Orbitals*, const DirectConfig&)
{
    ModuleBase::WARNING_QUIT("direct_gamma", "device_dm gpu requires real Gamma FP64 eigenvectors");
    return;
}

void direct_report_diagonal_start(const bool root, const std::chrono::steady_clock::time_point& started)
{
    const std::chrono::steady_clock::duration epoch = started.time_since_epoch();
    const long long start_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(epoch).count();
    if (root) std::printf("DM_DIRECT_DIAGONAL_START ns=%lld\n", start_ns);
}

void direct_report_redistribution(const bool root, const std::chrono::steady_clock::time_point& started)
{
    const double seconds = direct_elapsed(started);
    if (root) std::printf("DM_DIRECT_REDISTRIBUTION kind=hs seconds=%.9e\n", seconds);
}

}
