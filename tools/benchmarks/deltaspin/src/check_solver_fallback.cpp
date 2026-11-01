#include "module_hsolver/hsolver_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/kernels/dngvd_op.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_base/global_variable.h"
#include "module_base/parallel_global.h"
#include <base/macros/macros.h>
#include <cusolverDn.h>
#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;
using GPU = base_device::DEVICE_GPU;
using CPU = base_device::DEVICE_CPU;
using Solver = hsolver::HSolverPW<Complex, GPU>;

namespace
{
int current_point = 0, failed_point = -1, injections = 0;
int invalidations = 0, recovery_actions = 0;

class DiagonalOperator : public hamilt::Operator<Complex, GPU>
{
  public:
    DiagonalOperator() { this->cal_type = hamilt::calculation_type::pw_ekinetic; }

    void act(int columns, int leading, int polarizations, const Complex* input,
             Complex* output, int) const override
    {
        const size_t count = size_t(columns) * leading / polarizations;
        std::vector<Complex> values(count);
        cudaErrcheck(cudaMemcpy(values.data(), input, count * sizeof(Complex), cudaMemcpyDeviceToHost));
        for (size_t i = 0; i < count; ++i) { values[i] *= 1.0 + i % leading; }
        cudaErrcheck(cudaMemcpy(output, values.data(), count * sizeof(Complex), cudaMemcpyHostToDevice));
        if (invalidations != 0) { ++recovery_actions; }
    }
};

class DiagonalHamiltonian : public hamilt::Hamilt<Complex, GPU>
{
  public:
    void updateHk(int ik) override
    {
        current_point = ik;
        this->ops->init(ik);
    }
};

void run(int points, int polarizations)
{
    constexpr int plane_waves = 16, bands = 3;
    const int rows = plane_waves * polarizations;
    GlobalV::NPOL = polarizations;
    GlobalV::NSPIN = polarizations == 2 ? 4 : 1;
    ModulePW::PW_Basis_K basis;
    basis.nks = points;
    basis.npwk_max = plane_waves;
    basis.npwk = new int[points];
    std::fill(basis.npwk, basis.npwk + points, plane_waves);
    basis.gk2 = new double[points * plane_waves];
    std::fill(basis.gk2, basis.gk2 + points * plane_waves, 2.0);
    basis.tpiba2 = 1;
    K_Vectors kpoints;
    kpoints.wk.assign(points, 1.0 / points);
    elecstate::ElecState state;
    state.klist = &kpoints;
    state.ekb.create(points, bands);
    state.wg.create(points, bands);
    std::fill(state.wg.c, state.wg.c + points * bands, 1.0 / points);
    psi::Psi<Complex, GPU> psi(points, bands, rows, basis.npwk);
    std::vector<Complex> input(points * bands * rows);
    for (int ik = 0; ik < points; ++ik)
    for (int band = 0; band < bands; ++band)
    {
        input[(ik * bands + band) * rows + band] = std::sqrt(0.8);
        input[(ik * bands + band) * rows + band + bands] = Complex(0, std::sqrt(0.2));
    }
    DiagonalOperator op;
    DiagonalHamiltonian hamiltonian;
    hamiltonian.ops = &op;
    Solver solver(&basis, nullptr);
    Solver::BasisObserver observer;
    observer.capture = [](int, int, int, bool) {};
    observer.rotate = [](int, int, const Complex*, int, bool) {};
    observer.invalidate = [](int ik) {
        if (ik != failed_point) { std::exit(2); }
        ++invalidations;
    };
    double maximum_error = 0;
    for (int call = 0; call < 3; ++call)
    {
        psi.fix_k(0);
        cudaErrcheck(cudaMemcpy(psi.get_pointer(), input.data(), input.size() * sizeof(Complex), cudaMemcpyHostToDevice));
        failed_point = call == 1 ? points - 1 : -1;
        injections = invalidations = recovery_actions = 0;
        solver.solve(&hamiltonian, psi, &state, "dav_subspace", true,
                     nullptr, nullptr, {}, &observer);
        if (call == 1 && (injections != 1 || invalidations != 1 || recovery_actions == 0))
        {
            std::cerr << "CG recovery missing: injections=" << injections
                      << " invalidations=" << invalidations << " actions=" << recovery_actions << '\n';
            std::exit(3);
        }
        if (call != 1 && (injections || invalidations)) { std::exit(4); }
        psi.fix_k(0);
        std::vector<Complex> result(input.size());
        cudaErrcheck(cudaMemcpy(result.data(), psi.get_pointer(), result.size() * sizeof(Complex), cudaMemcpyDeviceToHost));
        for (int ik = 0; ik < points; ++ik)
        for (int band = 0; band < bands; ++band)
        {
            double norm = 0, residual = 0;
            for (int row = 0; row < rows; ++row)
            {
                const auto value = result[(ik * bands + band) * rows + row];
                norm += std::norm(value);
                residual += std::norm((row + 1.0 - state.ekb(ik, band)) * value);
            }
            const double error = std::max({std::abs(state.ekb(ik, band) - band - 1),
                                          std::abs(norm - 1), std::sqrt(residual)});
            if (!std::isfinite(error) || error > 1e-7) { std::exit(5); }
            maximum_error = std::max(maximum_error, error);
        }
        std::cout << "points=" << points << " npol=" << polarizations << " call=" << call
                  << " injected=" << injections << " invalidated=" << invalidations
                  << " recovery_actions=" << recovery_actions << " eigenpair_error=" << maximum_error << '\n';
    }
}
} // namespace

extern "C" cusolverStatus_t __real_cusolverDnZhegvd(cusolverDnHandle_t, cusolverEigType_t,
    cusolverEigMode_t, cublasFillMode_t, int, cuDoubleComplex*, int, cuDoubleComplex*, int,
    double*, cuDoubleComplex*, int, int*);

extern "C" cusolverStatus_t __wrap_cusolverDnZhegvd(cusolverDnHandle_t handle, cusolverEigType_t type,
    cusolverEigMode_t mode, cublasFillMode_t uplo, int n, cuDoubleComplex* a, int lda,
    cuDoubleComplex* b, int ldb, double* values, cuDoubleComplex* work, int length, int* info)
{
    const auto status = __real_cusolverDnZhegvd(handle, type, mode, uplo, n, a, lda,
                                              b, ldb, values, work, length, info);
    if (status == CUSOLVER_STATUS_SUCCESS && current_point == failed_point && injections == 0)
    {
        const int failure = 1;
        cudaErrcheck(cudaMemcpy(info, &failure, sizeof(int), cudaMemcpyHostToDevice));
        ++injections;
    }
    return status;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    POOL_WORLD = MPI_COMM_WORLD;
    GlobalV::BASIS_TYPE = "pw";
    GlobalV::CALCULATION = "nscf";
    GlobalV::psi_initializer = true;
    GlobalV::PW_DIAG_NDIM = 2;
    hsolver::DiagoIterAssist<Complex, GPU>::PW_DIAG_THR = 1e-12;
    hsolver::DiagoIterAssist<Complex, GPU>::PW_DIAG_NMAX = 100;
    hsolver::DiagoIterAssist<Complex, GPU>::need_subspace = false;
    hsolver::createGpuBlasHandle();
    hsolver::createGpuSolverHandle();
    for (int points : {1, 2})
    for (int polarizations : {1, 2}) { run(points, polarizations); }
    hsolver::destroyGpuSolverHandle();
    hsolver::destoryBLAShandle();
    MPI_Finalize();
}
