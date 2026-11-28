#include "dm_cusolver.h"
#include "kernels/cuda/diag_cusolver.cuh"
#include "source_base/module_device/device.h"
#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/timer.h"
#include "source_base/tool_quit.h"
#include "source_estate/elecstate_tools.h"
#include <algorithm>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

namespace hsolver
{
namespace
{
bool gpu_dm_supported(const int context, const int kpoints, const CusolverDmConfig& config)
{
    if (config.method != "cusolver" || config.kpar != 1 || config.nspin != 1
        || kpoints != 1 || config.skip_charge)
    {
        return false;
    }
    const base_device::DeviceContext& device_context = base_device::DeviceContext::instance();
    int unsupported_rank = device_context.get_device_count() != 1
                           || device_context.get_local_rank() != config.world_rank;
    char scope = 'A';
    char topology = ' ';
    Cigsum2d(context, &scope, &topology, 1, 1, &unsupported_rank, 1, -1, -1);
    return unsupported_rank == 0;
}

void dm_from_cusolver(Diag_Cusolver_gvd* solver, const double* weights, const int bands,
                     double* host_dm, double* distributed_dm, const Parallel_Blkcyclic& layout,
                     int* full_desc)
{
    ModuleBase::timer::start("HSolverLCAO", "dm_from_cusolver");
    const bool root = layout.irproc() == 0 && layout.icproc() == 0;
    const bool serial = layout.nrproc() == 1 && layout.ncproc() == 1;
    const int context = layout.ictxt();
    int success = 1;
    std::string error;
    if (root)
    {
        success = solver->density_from_retained(bands, weights, host_dm, error);
    }
    if (!serial)
    {
        char scope = 'A';
        char topology = ' ';
        if (root) Cigebs2d(context, &scope, &topology, 1, 1, &success, 1);
        else Cigebr2d(context, &scope, &topology, 1, 1, &success, 1, 0, 0);
    }
    if (!success)
    {
        const MPI_Comm comm = layout.comm();
        const Parallel::ParaWorld world = Parallel::ParaWorld::make_mpi("direct_dm", comm);
        Parallel::bcast_string(error, world, 0);
        ModuleBase::WARNING_QUIT("dm_from_cusolver", error, EXIT_FAILURE);
    }
    if (!serial)
    {
        const int n = layout.nrglb();
        int distributed_desc[9];
        const int* descriptor = layout.desc();
        const int* descriptor_end = descriptor + 9;
        std::copy(descriptor, descriptor_end, distributed_desc);
        // Bound ScaLAPACK's temporary packing buffers to one column block.
        const int block_size = layout.nb();
        for (int column = 0; column < n; column += block_size)
        {
            const int columns = std::min(block_size, n - column);
            const int first = column + 1;
            Cpdgemr2d(n, columns, host_dm, 1, first, full_desc,
                      distributed_dm, 1, first, distributed_desc, context);
        }
    }
    ModuleBase::timer::end("HSolverLCAO", "dm_from_cusolver");
}
}

bool solve_with_gpu_dm(HSMatrix<double>& hs, psi::Psi<double>& wave, elecstate::ElecState* state,
                  module_dm::DensityMatrix<double, double>& density, const Parallel_Orbitals* orbitals,
                  const CusolverDmConfig& config)
{
    ModuleBase::timer::start("HSolverLCAO", "solve_with_gpu_dm");
    const int bands = config.bands;
    const Parallel_Blkcyclic& layout = orbitals->para_2d();
    const int context = layout.ictxt();
    const int kpoints = wave.get_nk();
    const bool supported = gpu_dm_supported(context, kpoints, config);
    if (!supported)
    {
        ModuleBase::timer::end("HSolverLCAO", "solve_with_gpu_dm");
        return false;
    }
    const int n = layout.nrglb();
    const bool root = layout.irproc() == 0 && layout.icproc() == 0;
    const bool serial = layout.nrproc() == 1 && layout.ncproc() == 1;
    const std::size_t matrix_size = static_cast<std::size_t>(n) * n;
    const int block_size = layout.nb();
    Parallel_Blkcyclic gathered_layout;
    int gathered_context = -1;
    if (!serial && root)
    {
        // A one-rank grid retains the original block size during redistribution.
        gathered_layout.init(n, n, block_size, MPI_COMM_SELF, false);
        gathered_context = gathered_layout.ictxt();
    }
    int full_desc[9] = {1, gathered_context, n, n, block_size, block_size, 0, 0, n};
    int wave_desc[9] = {1, gathered_context, n, bands, block_size, block_size, 0, 0, n};
    int distributed_desc[9];
    int distributed_wave_desc[9];
    const int* matrix_descriptor = layout.desc();
    const int* matrix_descriptor_end = matrix_descriptor + 9;
    const int* wave_descriptor = orbitals->desc_wfc;
    const int* wave_descriptor_end = wave_descriptor + 9;
    std::copy(matrix_descriptor, matrix_descriptor_end, distributed_desc);
    std::copy(wave_descriptor, wave_descriptor_end, distributed_wave_desc);
    ModuleBase::MatrixBlock<double> hk;
    ModuleBase::MatrixBlock<double> sk;
    hs.hs_at_k(0, hk, sk);
    wave.fix_k(0);
    double* wave_pointer = wave.get_pointer();
    double* eigenvalues = &state->ekb(0, 0);
    double* dm_pointer = density.get_dmk_ptr(0);
    double* host_wave_pointer = wave_pointer;
    double* host_dm = dm_pointer;
    const double* h_pointer = hk.p;
    const double* s_pointer = sk.p;
    std::vector<double> h;
    std::vector<double> s;
    if (!serial)
    {
        const std::size_t root_matrix_size = root ? matrix_size : 1;
        h.resize(root_matrix_size);
        s.resize(root_matrix_size);
        Cpdgemr2d(n, n, hk.p, 1, 1, distributed_desc, h.data(), 1, 1, full_desc, context);
        Cpdgemr2d(n, n, sk.p, 1, 1, distributed_desc, s.data(), 1, 1, full_desc, context);
        h_pointer = h.data();
        s_pointer = s.data();
        // After upload, H's gather buffer can receive eigenvectors and then DM.
        host_wave_pointer = h.data();
        host_dm = h.data();
    }
    std::unique_ptr<Diag_Cusolver_gvd> solver;
    if (root)
    {
        solver.reset(new Diag_Cusolver_gvd());
        bool solved = false;
        const bool subset = Diag_Cusolver_gvd::use_subset(n, bands);
        if (subset)
        {
            solved = solver->solve_subset(n, bands, h_pointer, s_pointer, eigenvalues, host_wave_pointer, true);
        }
        if (!solved)
        {
            solver->solve_full(n, bands, h_pointer, s_pointer, eigenvalues, host_wave_pointer, true);
        }
    }
    if (!serial)
    {
        std::vector<double>().swap(s);
        char scope = 'A';
        char topology = ' ';
        if (root) Cdgebs2d(context, &scope, &topology, bands, 1, eigenvalues, bands);
        else Cdgebr2d(context, &scope, &topology, bands, 1, eigenvalues, bands, 0, 0);
        Cpdgemr2d(n, bands, host_wave_pointer, 1, 1, wave_desc,
                  wave_pointer, 1, 1, distributed_wave_desc, context);
    }
    // Keep the owning solver alive through ABACUS's actual occupation calculation.
    elecstate::calculate_weights(state->ekb, state->wg, state->klist, state->eferm,
                                 state->f_en, state->nelec_spin, bands, state->skip_weights);
    elecstate::calEBand(state->ekb, state->wg, state->f_en);
    Diag_Cusolver_gvd* retained = solver.get();
    const double* weights = &state->wg(0, 0);
    dm_from_cusolver(retained, weights, bands, host_dm, dm_pointer, layout, full_desc);
    ModuleBase::timer::end("HSolverLCAO", "solve_with_gpu_dm");
    return true;
}

bool solve_with_gpu_dm(HSMatrix<std::complex<double>>&, psi::Psi<std::complex<double>>&,
                  elecstate::ElecState*, module_dm::DensityMatrix<std::complex<double>, double>&,
                  const Parallel_Orbitals*, const CusolverDmConfig&)
{
    return false;
}

}
