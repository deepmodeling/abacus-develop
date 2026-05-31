/**
 * @file bench_real_comm.cpp
 * @brief Benchmark using ACTUAL ABACUS gather/scatter code from both branches.
 *
 * This benchmark directly calls PW_Basis::gatherp_scatters() / gathers_scatterp()
 * (feat/unblock nonblocking) and compares against the exact blocking implementations
 * extracted from the develop branch.
 *
 * Compile from abacus-develop root:
 *   g++ -std=c++14 -O3 -fopenmp -D__MPI \
 *     -I./source -I./source/source_basis/module_pw \
 *     -I./source/source_base -I./source/source_base/module_container \
 *     -I/usr/lib/x86_64-linux-gnu/openmpi/include \
 *     -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi \
 *     -o bench_real_comm \
 *     source/source_basis/module_pw/test/bench_real_comm.cpp \
 *     source/source_basis/module_pw/pw_basis.cpp \
 *     source/source_basis/module_pw/pw_init.cpp \
 *     source/source_basis/module_pw/pw_distributeg.cpp \
 *     source/source_basis/module_pw/pw_distributeg_method1.cpp \
 *     source/source_basis/module_pw/pw_distributeg_method2.cpp \
 *     source/source_basis/module_pw/pw_distributer.cpp \
 *     source/source_basis/module_pw/pw_basis_k.cpp \
 *     source/source_basis/module_pw/pw_basis_sup.cpp \
 *     source/source_basis/module_pw/pw_transform.cpp \
 *     source/source_basis/module_pw/pw_transform_k.cpp \
 *     source/source_basis/module_pw/pw_transform_gpu.cpp \
 *     source/source_base/module_fft/fft_bundle.cpp \
 *     source/source_base/module_fft/fft_cpu.cpp \
 *     source/source_base/mymath.cpp \
 *     source/source_base/timer.cpp \
 *     source/source_base/memory_recorder.cpp \
 *     source/source_base/tool_quit.cpp \
 *     source/source_base/matrix.cpp \
 *     source/source_base/matrix3.cpp \
 *     source/source_base/complexmatrix.cpp \
 *     source/source_base/module_external/blas_connector_base.cpp \
 *     source/source_base/module_external/blas_connector_vector.cpp \
 *     source/source_base/module_external/blas_connector_matrix.cpp \
 *     source/source_base/libm/branred.cpp \
 *     source/source_base/libm/sincos.cpp \
 *     source/source_base/module_device/memory_op.cpp \
 *     source/source_base/global_variable.cpp \
 *     -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi -lmpi_cxx \
 *     -lfftw3 -lopenblas -llapack -lpthread
 *
 * Run:
 *   OMP_NUM_THREADS=2 mpirun -np 4 ./bench_real_comm 500 64 64 64
 */

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <vector>

#include "../pw_basis.h"

// Global variables needed by ABACUS test infrastructure
int nproc_in_pool = 1;
int rank_in_pool = 0;
std::string precision_flag = "double";
std::string device_flag = "cpu";

// =========================================================================
//  EXACT blocking implementations from develop branch
//  (extracted from git show develop:source/source_basis/module_pw/pw_gatherscatter.h)
// =========================================================================

namespace BlockingOriginal {

template <typename T>
void gatherp_scatters(const ModulePW::PW_Basis& pw,
                      std::complex<T>* in,
                      std::complex<T>* out)
{
    if (pw.poolnproc == 1) {
        const int nst_ = pw.nst;
        const int nz_ = pw.nz;
        const int* istot2ixy_ = pw.istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int is = 0; is < nst_; ++is) {
            int ixy = istot2ixy_[is];
            std::complex<T>* outp = &out[is * nz_];
            std::complex<T>* inp  = &in[ixy * nz_];
            for (int iz = 0; iz < nz_; ++iz) {
                outp[iz] = inp[iz];
            }
        }
        return;
    }

#ifdef __MPI
    // Pack: (nplane fftnxy) to (nplane, nstot)
    const int nstot_gps = pw.nstot;
    const int nplane_gps = pw.nplane;
    const int* istot2ixy_gps = pw.istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int istot = 0; istot < nstot_gps; ++istot) {
        int ixy = istot2ixy_gps[istot];
        std::complex<T>* outp = &out[istot * nplane_gps];
        std::complex<T>* inp  = &in[ixy * nplane_gps];
        for (int iz = 0; iz < nplane_gps; ++iz) {
            outp[iz] = inp[iz];
        }
    }

    // Blocking Alltoallv
    if (typeid(T) == typeid(double)) {
        MPI_Alltoallv(out, pw.numr, pw.startr, MPI_DOUBLE_COMPLEX,
                      in,  pw.numg, pw.startg, MPI_DOUBLE_COMPLEX,
                      pw.pool_world);
    } else if (typeid(T) == typeid(float)) {
        MPI_Alltoallv(out, pw.numr, pw.startr, MPI_COMPLEX,
                      in,  pw.numg, pw.startg, MPI_COMPLEX,
                      pw.pool_world);
    }

    // Unpack: (numz[ip], ns) to (nz, ns)
    const int poolnproc_gps = pw.poolnproc;
    const int nst_gps = pw.nst;
    const int nz_gps = pw.nz;
    const int* numz_gps = pw.numz;
    const int* startg_gps = pw.startg;
    const int* startz_gps = pw.startz;
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < poolnproc_gps; ++ip) {
        for (int is = 0; is < nst_gps; ++is) {
            int nzip = numz_gps[ip];
            std::complex<T>* outp0 = &out[startz_gps[ip]];
            std::complex<T>* inp0  = &in[startg_gps[ip]];
            std::complex<T>* outp  = &outp0[is * nz_gps];
            std::complex<T>* inp   = &inp0[is * nzip];
            for (int izip = 0; izip < nzip; ++izip) {
                outp[izip] = inp[izip];
            }
        }
    }
#endif
}

template <typename T>
void gathers_scatterp(const ModulePW::PW_Basis& pw,
                      std::complex<T>* in,
                      std::complex<T>* out)
{
    if (pw.poolnproc == 1) {
        const int nrxx_ = pw.nrxx;
        const int nst_ = pw.nst;
        const int nz_ = pw.nz;
        const int* istot2ixy_ = pw.istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < nrxx_; ++i) {
            out[i] = std::complex<T>(0, 0);
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int is = 0; is < nst_; ++is) {
            int ixy = istot2ixy_[is];
            std::complex<T>* outp = &out[ixy * nz_];
            std::complex<T>* inp  = &in[is * nz_];
            for (int iz = 0; iz < nz_; ++iz) {
                outp[iz] = inp[iz];
            }
        }
        return;
    }

#ifdef __MPI
    // Pack: (nz, ns) to (numz[ip], ns, poolnproc)
    const int poolnproc_ = pw.poolnproc;
    const int nst_ = pw.nst;
    const int nz_ = pw.nz;
    const int* numz_ = pw.numz;
    const int* startg_ = pw.startg;
    const int* startz_ = pw.startz;
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int ip = 0; ip < poolnproc_; ++ip) {
        for (int is = 0; is < nst_; ++is) {
            int nzip = numz_[ip];
            std::complex<T>* outp0 = &out[startg_[ip]];
            std::complex<T>* inp0  = &in[startz_[ip]];
            std::complex<T>* outp  = &outp0[is * nzip];
            std::complex<T>* inp   = &inp0[is * nz_];
            for (int izip = 0; izip < nzip; ++izip) {
                outp[izip] = inp[izip];
            }
        }
    }

    // Blocking Alltoallv
    if (typeid(T) == typeid(double)) {
        MPI_Alltoallv(out, pw.numg, pw.startg, MPI_DOUBLE_COMPLEX,
                      in,  pw.numr, pw.startr, MPI_DOUBLE_COMPLEX,
                      pw.pool_world);
    } else if (typeid(T) == typeid(float)) {
        MPI_Alltoallv(out, pw.numg, pw.startg, MPI_COMPLEX,
                      in,  pw.numr, pw.startr, MPI_COMPLEX,
                      pw.pool_world);
    }

    // Unpack: (numz[ip], ns) to (nplane, nstot) and then to (nplane fftnxy)
    const int nrxx_gsp = pw.nrxx;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < nrxx_gsp; ++i) {
        out[i] = std::complex<T>(0, 0);
    }

    const int nstot = pw.nstot;
    const int nplane = pw.nplane;
    const int* istot2ixy = pw.istot2ixy;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int istot = 0; istot < nstot; ++istot) {
        int ixy = istot2ixy[istot];
        std::complex<T>* outp = &out[ixy * nplane];
        std::complex<T>* inp  = &in[istot * nplane];
        for (int iz = 0; iz < nplane; ++iz) {
            outp[iz] = inp[iz];
        }
    }
#endif
}

} // namespace BlockingOriginal

// =========================================================================
//  Benchmark using REAL ABACUS PW_Basis
// =========================================================================

class PWBasisAccessor : public ModulePW::PW_Basis {
public:
    PWBasisAccessor(const std::string& device_, const std::string& precision_)
        : ModulePW::PW_Basis(device_, precision_) {}

    using ModulePW::PW_Basis::gatherp_scatters;
    using ModulePW::PW_Basis::gathers_scatterp;
    using ModulePW::PW_Basis::distribute_r;
    using ModulePW::PW_Basis::distribute_g;
    using ModulePW::PW_Basis::getstartgr;
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n_iter = 500;
    if (argc > 1) n_iter = atoi(argv[1]);

    // Setup PW_Basis with realistic ABACUS parameters.
    // Use a 10-Angstrom cubic cell with ecut=100 Ry (typical for production runs).
    PWBasisAccessor pwtest("cpu", "double");
    const double lat0 = 1.8897261254578281;  // Bohr per Angstrom
    ModuleBase::Matrix3 latvec(10.0, 0.0, 0.0,
                                0.0, 10.0, 0.0,
                                0.0, 0.0, 10.0);
    const double gridecut = 100.0;  // Ry, determines FFT grid
    const double wfcecut  = 100.0;  // Ry, plane wave cutoff

#ifdef __MPI
    pwtest.initmpi(size, rank, MPI_COMM_WORLD);
#endif
    // Let ABACUS determine grid from ecut (more realistic), matching real production runs
    pwtest.initgrids(lat0, latvec, gridecut);
    pwtest.initparameters(false, wfcecut, 2 /*distribution_type=2*/, true);
    // Call setup steps individually, skipping fft_bundle.setupFFT()
    // (FFT plans are not needed for gather/scatter benchmarking)
    pwtest.distribute_r();
    pwtest.distribute_g();
    pwtest.getstartgr();

    // Allocate buffers matching ABACUS sizes
    int nrxx = pwtest.nrxx;
    int nst = pwtest.nst;
    int nz_local = pwtest.nz;
    int nplane = pwtest.nplane;
    int nstot = pwtest.nstot;

    int64_t plane_sz = (nrxx > 0) ? nrxx : 1;
    int64_t sticks_sz = (nst * nz_local > 0) ? nst * nz_local : 1;

    std::vector<std::complex<double>> plane_in(plane_sz);
    std::vector<std::complex<double>> sticks(sticks_sz);
    std::vector<std::complex<double>> plane_out(plane_sz);

    // Fill deterministic test data
    for (int64_t i = 0; i < plane_sz; ++i)
        plane_in[i] = std::complex<double>(sin(i * 0.01 + rank), cos(i * 0.01 + rank));

    if (rank == 0) {
        printf("=== REAL ABACUS gather/scatter benchmark ===\n");
        printf("Grid: %dx%dx%d (auto)  MPI ranks: %d  OMP threads: %d  Iterations: %d\n",
               pwtest.nx, pwtest.ny, pwtest.nz, size, omp_get_max_threads(), n_iter);
        printf("PW_Basis: nst=%d nz=%d nplane=%d nrxx=%d nstot=%d\n",
               nst, nz_local, nplane, nrxx, nstot);
        printf("numg_total=%d numr_total=%d\n",
               pwtest.startg[size-1] + pwtest.numg[size-1],
               pwtest.startr[size-1] + pwtest.numr[size-1]);
    }

    int warmup = std::max(10, n_iter / 10);

    // ---- Warmup nonblocking (feat/unblock, using actual PW_Basis) ----
    for (int i = 0; i < warmup; ++i) {
        pwtest.gatherp_scatters(plane_in.data(), sticks.data());
        pwtest.gathers_scatterp(sticks.data(), plane_out.data());
    }

    // ---- Warmup blocking (develop branch, exact code) ----
    for (int i = 0; i < warmup; ++i) {
        BlockingOriginal::gatherp_scatters<double>(pwtest, plane_in.data(), sticks.data());
        BlockingOriginal::gathers_scatterp<double>(pwtest, sticks.data(), plane_out.data());
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ---- TIMED: Nonblocking gatherp (feat/unblock) ----
    double t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        pwtest.gatherp_scatters(plane_in.data(), sticks.data());
    double t_nb_fwd = MPI_Wtime() - t0;

    // ---- TIMED: Blocking gatherp (develop) ----
    t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        BlockingOriginal::gatherp_scatters<double>(pwtest, plane_in.data(), sticks.data());
    double t_blk_fwd = MPI_Wtime() - t0;

    // ---- TIMED: Nonblocking gathers (feat/unblock) ----
    t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        pwtest.gathers_scatterp(sticks.data(), plane_out.data());
    double t_nb_rev = MPI_Wtime() - t0;

    // ---- TIMED: Blocking gathers (develop) ----
    t0 = MPI_Wtime();
    for (int i = 0; i < n_iter; ++i)
        BlockingOriginal::gathers_scatterp<double>(pwtest, sticks.data(), plane_out.data());
    double t_blk_rev = MPI_Wtime() - t0;

    // Collect stats
    double local[8], global_sum[8], global_max[8];
    local[0] = t_blk_fwd; local[1] = t_nb_fwd;
    local[2] = t_blk_rev; local[3] = t_nb_rev;
    local[4] = t_blk_fwd; local[5] = t_nb_fwd;
    local[6] = t_blk_rev; local[7] = t_nb_rev;

    MPI_Reduce(local, global_sum, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local + 4, global_max, 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double fwd_blk = global_sum[0] / size / n_iter;
        double fwd_nb  = global_sum[1] / size / n_iter;
        double rev_blk = global_sum[2] / size / n_iter;
        double rev_nb  = global_sum[3] / size / n_iter;

        printf("\n--- Per-iteration (us), average across ranks ---\n");
        printf("%-30s %10s %10s %8s\n", "Operation", "Blocking", "Nonblock", "Speedup");
        printf("%-30s %10.2f %10.2f %7.2fx\n",
               "gatherp_scatters (orig impl)", fwd_blk * 1e6, fwd_nb * 1e6, fwd_blk / fwd_nb);
        printf("%-30s %10.2f %10.2f %7.2fx\n",
               "gathers_scatterp (orig impl)", rev_blk * 1e6, rev_nb * 1e6, rev_blk / rev_nb);

        double tot_blk = fwd_blk + rev_blk;
        double tot_nb  = fwd_nb + rev_nb;
        printf("%-30s %10.2f %10.2f %7.2fx\n",
               "TOTAL roundtrip", tot_blk * 1e6, tot_nb * 1e6, tot_blk / tot_nb);

        printf("\n--- Max (slowest rank) per-iteration (us) ---\n");
        printf("gatherp  blocking=%10.2f  nonblocking=%10.2f  speedup=%.2fx\n",
               global_max[0] / n_iter * 1e6, global_max[1] / n_iter * 1e6,
               global_max[0] / global_max[1]);
        printf("gathers  blocking=%10.2f  nonblocking=%10.2f  speedup=%.2fx\n",
               global_max[2] / n_iter * 1e6, global_max[3] / n_iter * 1e6,
               global_max[2] / global_max[3]);

        double max_tot_blk = (global_max[0] + global_max[2]) / n_iter;
        double max_tot_nb  = (global_max[1] + global_max[3]) / n_iter;
        printf("TOTAL     blocking=%10.2f  nonblocking=%10.2f  speedup=%.2fx\n",
               max_tot_blk * 1e6, max_tot_nb * 1e6, max_tot_blk / max_tot_nb);

        printf("\n=== Speedup: %.2fx ===\n", tot_blk / tot_nb);
        if (tot_blk / tot_nb > 1.02)
            printf("Nonblocking (feat/unblock) is FASTER\n");
        else if (tot_blk / tot_nb < 0.98)
            printf("Blocking (develop) is FASTER\n");
        else
            printf("Performance is comparable\n");
    }

    MPI_Finalize();
    return 0;
}
