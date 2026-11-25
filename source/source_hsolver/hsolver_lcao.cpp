#include "hsolver_lcao.h"

#include "source_base/matrix_block.h"

#ifdef __MPI
#include "diago_scalapack.h"
#include "source_base/module_external/scalapack_connector.h"
#endif

#include "diago_lapack.h"

#ifdef __CUSOLVERMP
#include "diago_cusolvermp.h"
#endif

#ifdef __ELPA
#include "diago_elpa.h"
#include "diago_elpa_native.h"
#endif

#ifdef __CUDA
#include "diago_cusolver.h"
#endif

#ifdef __PEXSI
#include "diago_pexsi.h"
#endif

#include "source_base/module_device/device.h"
#include "source_estate/elecstate_tools.h"
#include "source_base/memory_recorder.h"
#include "source_base/timer.h"
#include "source_estate/elecstate_lcao.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_hsolver/parallel_k2d.h"

#ifdef ABACUS_DM_DIRECT_PILOT
#include "source_base/module_external/blacs_connector.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#endif

namespace hsolver
{

#ifdef ABACUS_DM_DIRECT_PILOT
namespace
{
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
    return std::chrono::duration<double>(finished - started).count();
}

template <typename T>
bool direct_gamma(HSMatrix<T>&, psi::Psi<T>&, elecstate::ElecState*,
                  module_dm::DensityMatrix<T, double>&, const Parallel_Orbitals*,
                  const int, const bool, const char*)
{
    return false;
}

bool direct_gamma(HSMatrix<double>& hs, psi::Psi<double>& wave, elecstate::ElecState* state,
                  module_dm::DensityMatrix<double, double>& density, const Parallel_Orbitals* orbitals,
                  const int bands, const bool supported, const char* mode)
{
    if (!supported || mode == nullptr || std::strcmp(mode, "cpu") == 0)
    {
        return false;
    }
    const bool verify = std::strcmp(mode, "verify") == 0;
    if (!verify && std::strcmp(mode, "direct") != 0) return false;
    const Parallel_Blkcyclic& layout = orbitals->para_2d();
    const int n = layout.nrglb();
    const int context = layout.ictxt();
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
    std::vector<double> host_wave(root_wave_size);
    std::vector<double> eigenvalues(n);
    std::vector<double> host_dm(root_matrix_size);
    const int leading_dimension = layout.irproc() == 0 ? n : 1;
    int full_desc[9] = {1, context, n, n, n, n, 0, 0, leading_dimension};
    int wave_desc[9] = {1, context, n, bands, n, bands, 0, 0, leading_dimension};
    int distributed_desc[9];
    int distributed_wave_desc[9];
    std::copy(layout.desc(), layout.desc() + 9, distributed_desc);
    std::copy(orbitals->desc_wfc, orbitals->desc_wfc + 9, distributed_wave_desc);
    const std::chrono::steady_clock::time_point gather_start = std::chrono::steady_clock::now();
    ModuleBase::MatrixBlock<double> hk;
    ModuleBase::MatrixBlock<double> sk;
    hs.hs_at_k(0, hk, sk);
    Cpdgemr2d(n, n, hk.p, 1, 1, distributed_desc, h.data(), 1, 1, full_desc, context);
    Cpdgemr2d(n, n, sk.p, 1, 1, distributed_desc, s.data(), 1, 1, full_desc, context);
    direct_sync(context);
    const double gather_seconds = direct_elapsed(gather_start);
    const std::chrono::steady_clock::time_point diagonal_start = std::chrono::steady_clock::now();
    std::unique_ptr<Diag_Cusolver_gvd> solver;
    if (root)
    {
        solver.reset(new Diag_Cusolver_gvd());
        solver->Dngvd_double_retained(n, h.data(), s.data(), eigenvalues.data(), host_wave.data(), bands);
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
    int success = 1;
    if (root)
    {
        const double* weights = &state->wg(0, 0);
        success = solver->density_from_retained(bands, weights, host_dm.data());
    }
    if (root) Cigebs2d(context, &scope, &topology, 1, 1, &success, 1);
    else Cigebr2d(context, &scope, &topology, 1, 1, &success, 1, 0, 0);
    const double device_dm_seconds = direct_elapsed(dm_start);
    const std::chrono::steady_clock::time_point scatter_start = std::chrono::steady_clock::now();
    double* dm_pointer = density.get_dmk_ptr(0);
    if (success)
    {
        Cpdgemr2d(n, n, host_dm.data(), 1, 1, full_desc,
                  dm_pointer, 1, 1, distributed_desc, context);
    }
    else
    {
        // Host psi remains valid even when the optional DM allocation fails.
        module_dm::dm_from_psi(orbitals, state->wg, wave, density);
        if (root) std::printf("DM_DIRECT_FALLBACK reason=device_dm_failure\n");
    }
    direct_sync(context);
    const double scatter_seconds = direct_elapsed(scatter_start);
    if (root)
    {
        std::printf("DM_DIRECT_PHASE mode=%s n=%d bands=%d gather_hs=%.9e diagonal=%.9e psi_redistribute=%.9e occupations=%.9e device_dm=%.9e dm_redistribute=%.9e hs_redistributions=2 hs_logical_bytes=%zu psi_redistributions=1 psi_logical_bytes=%zu dm_redistributions=%d dm_logical_bytes=%zu\n",
                    mode, n, bands, gather_seconds, diagonal_seconds, psi_seconds, occupation_seconds,
                    device_dm_seconds, scatter_seconds, hs_bytes, wave_bytes, success, matrix_bytes);
    }
    if (verify && success)
    {
        if (root)
        {
            double trace = 0.0;
            double electrons = 0.0;
            for (int band = 0; band < bands; ++band) electrons += state->wg(0, band);
            for (int col = 0; col < n; ++col)
            {
                for (int row = 0; row <= col; ++row)
                {
                    const std::size_t index = static_cast<std::size_t>(col) * n + row;
                    const double multiplicity = row == col ? 1.0 : 2.0;
                    trace += multiplicity * host_dm[index] * s[index];
                }
            }
            const double error = std::abs(trace - electrons);
            std::printf("DM_DIRECT_TRACE trace_ds=%.17e electrons=%.17e abs_error=%.17e\n", trace, electrons, error);
            if (!std::isfinite(trace) || error > 1.0e-7)
                ModuleBase::WARNING_QUIT("direct_gamma", "Tr(DS) failed the 1e-7 electron tolerance");
        }
        const std::size_t count = layout.size();
        std::vector<double> gpu(dm_pointer, dm_pointer + count);
        module_dm::dm_from_psi(orbitals, state->wg, wave, density);
        double max_error = 0.0;
        double error2 = 0.0;
        double reference2 = 0.0;
        for (std::size_t index = 0; index < count; ++index)
        {
            const double error = gpu[index] - dm_pointer[index];
            const double absolute = std::abs(error);
            max_error = std::max(max_error, absolute);
            error2 += error * error;
            reference2 += dm_pointer[index] * dm_pointer[index];
            if (!std::isfinite(error) || absolute > 1.0e-10)
                ModuleBase::WARNING_QUIT("direct_gamma", "GPU DM failed the 1e-10 element tolerance");
        }
        const int row = layout.irproc();
        const int col = layout.icproc();
        std::printf("DM_DIRECT_VERIFY row=%d col=%d count=%zu max_abs=%.17e error2=%.17e reference2=%.17e\n",
                    row, col, count, max_error, error2, reference2);
        std::copy(gpu.begin(), gpu.end(), dm_pointer);
    }
    return true;
}
}
#endif

template <typename TK>
void HSolverLCAO<TK>::solve(HSMatrix<TK>& hs,
                                   psi::Psi<TK>& psi,
								   elecstate::ElecState* pes,
								   module_dm::DensityMatrix<TK, double>& dm, // mohan add 2025-11-03
								   Charge &chr,
                                   const int nspin,
                                   const double omega,
                                   const bool skip_charge)
{
    ModuleBase::TITLE("HSolverLCAO", "solve");
    ModuleBase::timer::start("HSolverLCAO", "solve");

    if (this->method != "pexsi")
    {
#ifdef ABACUS_DM_DIRECT_PILOT
        const char* trial_mode = std::getenv("ABACUS_DM_DIRECT_PILOT");
        const Parallel_Blkcyclic& trial_layout = this->ParaV->para_2d();
        const int trial_context = trial_layout.ictxt();
        const bool trial_root = trial_layout.irproc() == 0 && trial_layout.icproc() == 0;
        direct_sync(trial_context);
        const std::chrono::steady_clock::time_point combined_start = std::chrono::steady_clock::now();
        const base_device::DeviceContext& device_context = base_device::DeviceContext::instance();
        int unsupported_rank = device_context.get_device_count() != 1
                               || device_context.get_local_rank() != this->world_rank;
        char scope = 'A';
        char topology = ' ';
        Cigsum2d(trial_context, &scope, &topology, 1, 1, &unsupported_rank, 1, -1, -1);
        const bool supported = this->method == "cusolver" && this->kpar_lcao == 1
                               && nspin == 1 && psi.get_nk() == 1 && !skip_charge
                               && unsupported_rank == 0 && this->world_nproc > 1;
        const bool direct_done = direct_gamma(hs, psi, pes, dm, this->ParaV, this->nbands, supported, trial_mode);
        if (!direct_done)
        {
            if (trial_root && trial_mode != nullptr && std::strcmp(trial_mode, "cpu") != 0)
                std::printf("DM_DIRECT_FALLBACK reason=unsupported_scope\n");
#endif
    #ifdef __MPI
    #ifdef __CUDA
        if (this->method == "cusolver" && this->world_nproc > 1)
        {
            this->parakSolve_cusolver(hs, psi, pes);
        }else 
    #endif
        if (this->kpar_lcao > 1
            && (this->method == "genelpa" || this->method == "elpa" || this->method == "scalapack_gvx" || this->method == "lapack"))
        {
            this->parakSolve(hs, psi, pes, this->kpar_lcao, nspin);
        } else
    #endif
        if (this->kpar_lcao == 1)
        {
            /// Loop over k points for solve Hamiltonian to eigenpairs(eigenvalues and eigenvectors).
            for (int ik = 0; ik < psi.get_nk(); ++ik)
            {
                /// H(k) and S(k) are all the eigensolvers need from the Hamiltonian
                ModuleBase::MatrixBlock<TK> hk, sk;
                hs.hs_at_k(ik, hk, sk);

                /// find psi pointer for each k point
                psi.fix_k(ik);

                /// solve eigenvector and eigenvalue for H(k)
                this->hamiltSolvePsiK(hk, sk, psi, &(pes->ekb(ik, 0)));
            }
        }
        else
        {
            ModuleBase::WARNING_QUIT("HSolverLCAO::solve",
                                     "This method and KPAR setting is not supported for lcao basis in ABACUS!");
        }

        elecstate::calculate_weights(pes->ekb,
                                     pes->wg,
                                     pes->klist,
                                     pes->eferm,
                                     pes->f_en,
                                     pes->nelec_spin,
                                     this->nbands,
                                     pes->skip_weights);

        elecstate::calEBand(pes->ekb, pes->wg, pes->f_en);
#ifdef ABACUS_DM_DIRECT_PILOT
        const std::chrono::steady_clock::time_point cpu_dm_start = std::chrono::steady_clock::now();
#endif
        module_dm::dm_from_psi(this->ParaV, pes->wg, psi, dm);
#ifdef ABACUS_DM_DIRECT_PILOT
        direct_sync(trial_context);
        const double cpu_dm_seconds = direct_elapsed(cpu_dm_start);
        if (trial_root) std::printf("DM_DIRECT_CPU_DM seconds=%.9e\n", cpu_dm_seconds);
        }
        direct_sync(trial_context);
        const double combined_seconds = direct_elapsed(combined_start);
        if (trial_root && trial_mode != nullptr)
            std::printf("DM_DIRECT_COMBINED mode=%s direct=%d seconds=%.9e\n", trial_mode, direct_done, combined_seconds);
#endif
        dm.cal_dmr(-1);

        if (!skip_charge)
        {
            // compute charge density from density matrix, mohan update 20251024
            // delegate to ElecStateLCAO to keep the source_lcao dependency out of
            // source_hsolver (mirrors the pexsi branch below and the PW psiToRho path)
            dynamic_cast<elecstate::ElecStateLCAO<TK>*>(pes)->dmToRho(dm.get_dmr_vec(), nspin, &chr, omega);
        }
        else
        {
            // used in nscf calculation
        }
    }
    else if (this->method == "pexsi")
    {
#ifdef __PEXSI // other purification methods should follow this routine
        DiagoPexsi<TK> pe(ParaV, nspin, this->nlocal, this->nelec, this->world_nproc);
        for (int ik = 0; ik < psi.get_nk(); ++ik)
        {
            psi.fix_k(ik);
            ModuleBase::MatrixBlock<TK> hk, sk;
            hs.hs_at_k(ik, hk, sk);
            // solve eigenvector and eigenvalue for H(k)
            pe.diag(hk, sk, psi, nullptr);
        }
        auto _pes = dynamic_cast<elecstate::ElecStateLCAO<TK>*>(pes);
        pes->f_en.eband = pe.totalFreeEnergy;
        // maybe eferm could be dealt with in the future
        _pes->dm2rho(pe.DM, pe.EDM, &dm, omega);
#endif
    }

    ModuleBase::timer::end("HSolverLCAO", "solve");
    return;
}

template <typename T>
void HSolverLCAO<T>::hamiltSolvePsiK(ModuleBase::MatrixBlock<T>& hk,
                                     ModuleBase::MatrixBlock<T>& sk,
                                     psi::Psi<T>& psi,
                                     double* eigenvalue)
{
    ModuleBase::TITLE("HSolverLCAO", "hamiltSolvePsiK");
    ModuleBase::timer::start("HSolverLCAO", "hamiltSolvePsiK");

    if (this->method == "scalapack_gvx")
    {
#ifdef __MPI
        DiagoScalapack<T> sa(this->nlocal, this->nbands);
        sa.diag(hk, sk, psi, eigenvalue);
#endif
    }
#ifdef __ELPA
    else if (this->method == "genelpa")
    {
        DiagoElpa<T> el(this->nlocal, this->nbands);
        el.diag(hk, sk, psi, eigenvalue);
    }
    else if (this->method == "elpa")
    {
        DiagoElpaNative<T> el(this->nlocal, this->nbands, this->use_gpu);
        el.diag(hk, sk, psi, eigenvalue);
    }
#endif
#ifdef __CUDA
    else if (this->method == "cusolver")
    {
        // Note: This branch will only be executed in the single-process case
        DiagoCusolver<T> cu(this->nlocal, this->nbands);
        cu.diag(hk, sk, psi, eigenvalue);
    }
#ifdef __CUSOLVERMP
    else if (this->method == "cusolvermp")
    {
        DiagoCusolverMP<T> cm(this->nlocal, this->nbands);
        cm.diag(hk, sk, psi, eigenvalue);
    }
#endif
#endif
    else if (this->method == "lapack") // only for single core
    {
        DiagoLapack<T> la(this->nlocal, this->nbands);
        la.diag(hk, sk, psi, eigenvalue);
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method is not supported for lcao basis in ABACUS!");
    }

    ModuleBase::timer::end("HSolverLCAO", "hamiltSolvePsiK");
}

template <typename T>
void HSolverLCAO<T>::parakSolve(HSMatrix<T>& hs,
                                        psi::Psi<T>& psi,
                                        elecstate::ElecState* pes,
                                        const int kpar,
                                        const int nspin)
{
#ifdef __MPI
    ModuleBase::timer::start("HSolverLCAO", "parakSolve");
    Parallel_K2D<T> k2d;
    k2d.set_kpar(kpar);
    int nbands = this->ParaV->get_nbands();
    int nks = psi.get_nk();
    int nrow = this->ParaV->para_2d().nrglb();
    int nb2d = this->ParaV->para_2d().nb();
    k2d.set_para_env(psi.get_nk(), nrow, nb2d, this->world_nproc, this->world_rank, nspin);
    /// set psi_pool
    const int zero = 0;
    int coord_col = k2d.get_p2D_pool()->icproc();
    const int pool_grid_cols = k2d.get_p2D_pool()->ncproc();
    int ncol_bands_pool = numroc_(&(nbands), &(nb2d), &coord_col, &zero, &pool_grid_cols);
    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < k2d.get_pKpoints()->get_max_nks_pool(); ++ik)
    {
        // if nks is not equal to the number of k points in the pool
        std::vector<int> ik_kpar;
        int ik_avail = 0;
        for (int i = 0; i < k2d.get_kpar(); i++)
        {
            if (ik + k2d.get_pKpoints()->startk_pool[i] < nks && ik < k2d.get_pKpoints()->nks_pool[i])
            {
                ik_avail++;
            }
        }
        if (ik_avail == 0)
        {
            ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "ik_avail is 0!");
        }
        else
        {
            ik_kpar.resize(ik_avail);
            for (int i = 0; i < ik_avail; i++)
            {
                ik_kpar[i] = ik + k2d.get_pKpoints()->startk_pool[i];
            }
        }
        k2d.distribute_hsk(hs, ik_kpar, nrow);
        /// global index of k point
        int ik_global = ik + k2d.get_pKpoints()->startk_pool[k2d.get_my_pool()];
        const int pool_local_rows = k2d.get_p2D_pool()->nrloc();
        auto psi_pool = psi::Psi<T>(1, ncol_bands_pool, pool_local_rows, pool_local_rows, true);
        ModuleBase::Memory::record("HSolverLCAO::psi_pool", nrow * ncol_bands_pool * sizeof(T));
        if (ik_global < psi.get_nk() && ik < k2d.get_pKpoints()->nks_pool[k2d.get_my_pool()])
        {
            /// local psi in pool
            psi_pool.fix_k(0);
            ModuleBase::MatrixBlock<T> hk_pool = ModuleBase::MatrixBlock<T>{k2d.hk_pool.data(),
                                                                            (size_t)k2d.get_p2D_pool()->nrloc(),
                                                                            (size_t)k2d.get_p2D_pool()->ncloc(),
                                                                            k2d.get_p2D_pool()->desc()};
            ModuleBase::MatrixBlock<T> sk_pool = ModuleBase::MatrixBlock<T>{k2d.sk_pool.data(),
                                                                            (size_t)k2d.get_p2D_pool()->nrloc(),
                                                                            (size_t)k2d.get_p2D_pool()->ncloc(),
                                                                            k2d.get_p2D_pool()->desc()};
            /// solve eigenvector and eigenvalue for H(k)
            if (this->method == "scalapack_gvx")
            {
                DiagoScalapack<T> sa(this->nlocal, this->nbands);
                sa.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
            else if (this->method == "lapack")
            {
                DiagoLapack<T> la(this->nlocal, this->nbands);
                la.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
#ifdef __ELPA
            else if (this->method == "genelpa")
            {
                DiagoElpa<T> el(this->nlocal, this->nbands);
                el.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
            else if (this->method == "elpa")
            {
                DiagoElpaNative<T> el(this->nlocal, this->nbands, this->use_gpu);
                el.diag_pool(hk_pool, sk_pool, psi_pool, &(pes->ekb(ik_global, 0)), k2d.POOL_WORLD_K2D);
            }
#endif
            else
            {
                ModuleBase::WARNING_QUIT("HSolverLCAO::solve",
                                         "This type of eigensolver for k-parallelism diagnolization is not supported!");
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ModuleBase::timer::start("HSolverLCAO", "collect_psi");
        for (int ipool = 0; ipool < ik_kpar.size(); ++ipool)
        {
            int source = k2d.get_pKpoints()->get_startpro_pool(ipool);
            MPI_Bcast(&(pes->ekb(ik_kpar[ipool], 0)), nbands, MPI_DOUBLE, source, MPI_COMM_WORLD);
            constexpr int DESCRIPTOR_ENTRY_COUNT = 9;
            int desc_pool[DESCRIPTOR_ENTRY_COUNT];
            const int* const pool_descriptor = k2d.get_p2D_pool()->desc();
            std::copy(pool_descriptor, pool_descriptor + DESCRIPTOR_ENTRY_COUNT, desc_pool);
            if (k2d.get_my_pool() != ipool)
            {
                desc_pool[1] = -1;
            }
            psi.fix_k(ik_kpar[ipool]);
            Cpxgemr2d(nrow,
                      nbands,
                      psi_pool.get_pointer(),
                      1,
                      1,
                      desc_pool,
                      psi.get_pointer(),
                      1,
                      1,
                      const_cast<int*>(k2d.get_p2D_global()->desc()),
                      k2d.get_p2D_global()->ictxt());
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ModuleBase::timer::end("HSolverLCAO", "collect_psi");
    }
    k2d.unset_para_env();
    ModuleBase::timer::end("HSolverLCAO", "parakSolve");
#endif
}

#if defined (__MPI) && defined (__CUDA)
template <typename T>
void HSolverLCAO<T>::parakSolve_cusolver(HSMatrix<T>& hs,
                                            psi::Psi<T>& psi,
                                            elecstate::ElecState* pes)
{
    ModuleBase::timer::start("HSolverLCAO", "parakSolve");
    // GPU device is already bound by DeviceContext::init() in read_input.cpp
    auto& dev_ctx = base_device::DeviceContext::instance();
    const int local_rank = dev_ctx.get_local_rank();
    const int device_count = dev_ctx.get_device_count();

    const int world_rank = this->world_rank;
    const int world_size = this->world_nproc;

    // Determine if this process is active (assigned a dedicated GPU)
    // We enforce 1 process per GPU by checking local_rank < device_count.
    bool is_active = (local_rank < device_count);
    int is_active_int = is_active ? 1 : 0;

    std::vector<int> all_local_ranks(world_size);
    std::vector<int> all_is_active(world_size);

    MPI_Allgather(&local_rank, 1, MPI_INT, 
                  all_local_ranks.data(), 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&is_active_int, 1, MPI_INT, 
                  all_is_active.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int total_active_ranks = 0;
    int max_local_rank = 0;
    for (int i = 0; i < world_size; ++i) {
        if (all_is_active[i]) {
            total_active_ranks++;
            if(all_local_ranks[i] > max_local_rank) max_local_rank = all_local_ranks[i];
        }
    }

    // active_ranks will store the global ranks of all active processes across all nodes
    // The order of global ranks stored here determines the order in which they will be assigned to k-points.
    // The k-points will be distributed among these ranks in a round-robin fashion.
    // The purpose of setting the order is to ensure load balancing among nodes as much as possible
    std::vector<int> active_ranks;
    for(int r = 0; r <= max_local_rank; r++)
    {
        for(int j = 0; j < world_size; j++)
        {
            if(all_is_active[j] && all_local_ranks[j] == r)
            {
                active_ranks.push_back(j);
            }
        }
    }
    
    const int nks = psi.get_nk();  // total number of k points
    const int nbands = this->ParaV->get_nbands();
    // Set the parallel storage scheme for the matrix and psi
    Parallel_Blkcyclic mat_para_global;    // store the info about how the origin matrix is distributed in parallel
    Parallel_Blkcyclic mat_para_local;     // store the info about how the matrix is distributed after collected from all processes
    Parallel_Blkcyclic psi_para_global;    // store the info about how the psi is distributed in parallel
    Parallel_Blkcyclic psi_para_local;     // store the info about how the psi is distributed before distributing to all processes

    MPI_Comm self_comm;  // the communicator that only contains the current process itself
    MPI_Comm_split(MPI_COMM_WORLD, world_rank, 0, &self_comm);
    int nrow = this->ParaV->para_2d().nrglb(); // number of rows in the global matrix
    int ncol = nrow;
    int nb2d = this->ParaV->para_2d().nb();      // block size for the 2D matrix distribution
    mat_para_global.init(nrow, ncol, nb2d, MPI_COMM_WORLD, false);
    psi_para_global.init(nrow, nbands, nb2d, MPI_COMM_WORLD, false);
    mat_para_local.init(nrow, ncol, nb2d, self_comm, false);
    psi_para_local.init(nrow, ncol, nb2d, self_comm, false);

    std::vector<T> hk_mat; // temporary storage for H(k) matrix collected from all processes
    std::vector<T> sk_mat; // temporary storage for S(k) matrix collected from all processes
    constexpr int DESCRIPTOR_ENTRY_COUNT = 9;
    constexpr int DESCRIPTOR_BLACS_CONTEXT_INDEX = 1;
    // In each iteration, we process total_active_ranks k-points.
    for(int ik_start = 0; ik_start < nks; ik_start += total_active_ranks)
    {
#ifdef ABACUS_DM_DIRECT_PILOT
        const std::chrono::steady_clock::time_point gather_start = std::chrono::steady_clock::now();
#endif
        int kpt_assigned = -1; // the k-point assigned to the current MPI process in this iteration
        // Compute and gather the hk and sk matrices distributed across different processes in parallel,
        // preparing for subsequent transfer to the GPU for computation.
        for(int ik = ik_start; ik < ik_start + total_active_ranks && ik < nks; ik++)
        {
            // `is_active` indicates whether this MPI process is assigned to compute the current k-point
            bool is_active = world_rank == active_ranks[ik % total_active_ranks];
            if (is_active)
            {
                kpt_assigned = ik;
                hk_mat.resize(nrow * ncol);
                sk_mat.resize(nrow * ncol);
            }
            ModuleBase::MatrixBlock<T> hk_2D, sk_2D;
            hs.hs_at_k(ik, hk_2D, sk_2D);
            int desc_tmp[DESCRIPTOR_ENTRY_COUNT];
            T* hk_local_ptr = hk_mat.data();
            T* sk_local_ptr = sk_mat.data();
            const int* const local_descriptor = mat_para_local.desc();
            std::copy(local_descriptor, local_descriptor + DESCRIPTOR_ENTRY_COUNT, desc_tmp);
            if( !is_active)
            {
                desc_tmp[DESCRIPTOR_BLACS_CONTEXT_INDEX] = -1;
            }

            Cpxgemr2d(nrow, ncol, hk_2D.p, 1, 1,
                      const_cast<int*>(mat_para_global.desc()),
                      hk_local_ptr, 1, 1, desc_tmp,
                      mat_para_global.ictxt());
            Cpxgemr2d(nrow, ncol, sk_2D.p, 1, 1,
                      const_cast<int*>(mat_para_global.desc()),
                      sk_local_ptr, 1, 1, desc_tmp,
                      mat_para_global.ictxt());
        }

        // diagonalize the Hamiltonian matrix using cusolver
#ifdef ABACUS_DM_DIRECT_PILOT
        direct_sync(mat_para_global.ictxt());
        const double gather_seconds = direct_elapsed(gather_start);
        const std::chrono::steady_clock::time_point diagonal_start = std::chrono::steady_clock::now();
#endif
        psi::Psi<T> psi_local{};
        if(kpt_assigned != -1)
        {
            psi_local.resize(1, ncol, nrow);
            DiagoCusolver<T> cu(this->nlocal, this->nbands);
            ModuleBase::MatrixBlock<T> hk_local = ModuleBase::MatrixBlock<T>{
                    hk_mat.data(), (size_t)nrow, (size_t)ncol,
                    mat_para_local.desc()};
            ModuleBase::MatrixBlock<T> sk_local = ModuleBase::MatrixBlock<T>{
                    sk_mat.data(), (size_t)nrow, (size_t)ncol,
                    mat_para_local.desc()};
            cu.diag(hk_local, sk_local, psi_local, &(pes->ekb(kpt_assigned, 0)));
        }

        // transfer the eigenvectors and eigenvalues to all processes
#ifdef ABACUS_DM_DIRECT_PILOT
        direct_sync(mat_para_global.ictxt());
        const double diagonal_seconds = direct_elapsed(diagonal_start);
        const std::chrono::steady_clock::time_point scatter_start = std::chrono::steady_clock::now();
#endif
        for(int ik = ik_start; ik < ik_start + total_active_ranks && ik < nks; ik++)
        {
            int root = active_ranks[ik % total_active_ranks];
            MPI_Bcast(&(pes->ekb(ik, 0)), nbands, MPI_DOUBLE, root, MPI_COMM_WORLD);
            int desc_pool[DESCRIPTOR_ENTRY_COUNT];
            const int* const psi_local_descriptor = psi_para_local.desc();
            std::copy(psi_local_descriptor,
                      psi_local_descriptor + DESCRIPTOR_ENTRY_COUNT,
                      desc_pool);
            T* psi_local_ptr = nullptr;
            if (world_rank != root)
            {
                desc_pool[1] = -1;
            }else
            {
                psi_local_ptr = psi_local.get_pointer();
            }
            psi.fix_k(ik);
            Cpxgemr2d(nrow,
                      nbands,
                      psi_local_ptr,
                      1,
                      1,
                      desc_pool,
                      psi.get_pointer(),
                      1,
                      1,
                      const_cast<int*>(psi_para_global.desc()),
                      psi_para_global.ictxt());
        }
#ifdef ABACUS_DM_DIRECT_PILOT
        direct_sync(mat_para_global.ictxt());
        const double psi_seconds = direct_elapsed(scatter_start);
        if (world_rank == 0)
        {
            const std::size_t matrix_bytes = sizeof(T) * static_cast<std::size_t>(nrow) * nrow;
            const std::size_t wave_bytes = sizeof(T) * static_cast<std::size_t>(nrow) * nbands;
            const std::size_t hs_bytes = 2 * matrix_bytes;
            std::printf("DM_DIRECT_PHASE mode=cpu n=%d bands=%d gather_hs=%.9e diagonal=%.9e psi_redistribute=%.9e hs_redistributions=2 hs_logical_bytes=%zu psi_redistributions=1 psi_logical_bytes=%zu dm_redistributions=0 dm_logical_bytes=0\n",
                        nrow, nbands, gather_seconds, diagonal_seconds, psi_seconds, hs_bytes, wave_bytes);
        }
#endif
    }

    MPI_Comm_free(&self_comm);
    ModuleBase::timer::end("HSolverLCAO", "parakSolve");
}
#endif


template class HSolverLCAO<double>;
template class HSolverLCAO<std::complex<double>>;

} // namespace hsolver
