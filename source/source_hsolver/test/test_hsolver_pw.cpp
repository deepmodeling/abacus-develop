#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>

#include "hsolver_pw_sup.h"
#include "hsolver_supplementary_mock.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/diago_ppcg.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/hsolver_lcaopw.h"
#include "source_hsolver/hsolver_pw.h"
#include "source_basis/module_pw/test/test_tool.h"
#include "source_base/parallel_reduce.h"

#include <algorithm>
#include <complex>
#include <cmath>
#include <limits>

/// H = S = identity: the simplest operator the solvers can be handed
template <typename T>
class IdentityHSOperator : public hsolver::HSOperator<T, base_device::DEVICE_CPU>
{
  public:
    void update_k(const int ik) override
    {
    }
    void hpsi(const T* x, T* hx, const int ld, const int nvec) const override
    {
        std::copy(x, x + static_cast<size_t>(ld) * nvec, hx);
    }
    void spsi(const T* x, T* sx, const int ld, const int nvec) const override
    {
        std::copy(x, x + static_cast<size_t>(ld) * nvec, sx);
    }
};

class DiagonalPPCGOperator : public IdentityHSOperator<std::complex<double>>
{
  public:
    void hpsi(const std::complex<double>* x,
              std::complex<double>* hx,
              const int ld,
              const int nvec) const override
    {
        block_sizes.push_back(nvec);
        for (int ib = 0; ib < nvec; ++ib)
        {
            for (int ig = 0; ig < ld; ++ig)
            {
                hx[ib * ld + ig] = double(ig + 1) * x[ib * ld + ig];
            }
        }
    }

    mutable std::vector<int> block_sizes;
};

// Mock implementations for the template functions causing linking errors
namespace ModulePW {
    // Mock implementation for recip_to_real
    template<typename FPTYPE, typename Device>
    void PW_Basis_K::recip_to_real(const Device* ctx,
                                  const std::complex<FPTYPE>* in,
                                  std::complex<FPTYPE>* out,
                                  const int ik,
                                  const bool add,
                                  const FPTYPE factor) const
    {
        // Simple mock implementation that does nothing
        // In a real test, you might want to implement behavior that simulates the actual function
    }

    // Mock implementation for real_to_recip
    template<typename FPTYPE, typename Device>
    void PW_Basis_K::real_to_recip(const Device* ctx,
                                  const std::complex<FPTYPE>* in,
                                  std::complex<FPTYPE>* out,
                                  const int ik,
                                  const bool add,
                                  const FPTYPE factor) const
    {
        // Simple mock implementation that does nothing
    }

    // Explicit template instantiations
    template void PW_Basis_K::recip_to_real<float, base_device::DEVICE_CPU>(
        const base_device::DEVICE_CPU* ctx,
        const std::complex<float>* in,
        std::complex<float>* out,
        const int ik,
        const bool add,
        const float factor) const;

    template void PW_Basis_K::recip_to_real<double, base_device::DEVICE_CPU>(
        const base_device::DEVICE_CPU* ctx,
        const std::complex<double>* in,
        std::complex<double>* out,
        const int ik,
        const bool add,
        const double factor) const;

    template void PW_Basis_K::real_to_recip<float, base_device::DEVICE_CPU>(
        const base_device::DEVICE_CPU* ctx,
        const std::complex<float>* in,
        std::complex<float>* out,
        const int ik,
        const bool add,
        const float factor) const;

    template void PW_Basis_K::real_to_recip<double, base_device::DEVICE_CPU>(
        const base_device::DEVICE_CPU* ctx,
        const std::complex<double>* in,
        std::complex<double>* out,
        const int ik,
        const bool add,
        const double factor) const;
}

/************************************************
 *  unit test of HSolverPW class
 ***********************************************/

/**
 * Tested function:
 *  - test for template float and double respectively:
 *  - 1. solve()
 *  - 2. initDiagh()
 *  - 3. endDiagh()
 *  - 4. hamiltSolvePsiK()
 *  - 5. updatePsiK()
 *  - 6. update_precondition()
 *  - 7. hsolver::HSolver::diagethr (for cases below)
 * 		- set_diagethr, for setting diagethr;
 *  - 8. solve()
 *      - lcao_in_pw specific implementation
 */

// mock diago_hs_para
namespace hsolver {
template <typename T>
void diago_hs_para(T* h,
                   T* s,
                   const int lda,
                   const int nband,
                   typename GetTypeReal<T>::type* const ekb,
                   T* const wfc,
                   const MPI_Comm& comm,
                   const int diag_subspace,
                   const int block_size = 0)
{}
template void diago_hs_para<double>(double* h,
                                    double* s,
                                    const int lda,
                                    const int nband,
                                    typename GetTypeReal<double>::type* const ekb,
                                    double* const wfc,
                                    const MPI_Comm& comm,
                                    const int diag_subspace,
                                    const int block_size);

template void diago_hs_para<std::complex<double>>(std::complex<double>* h,
                                                  std::complex<double>* s,
                                                  const int lda,
                                                  const int nband,
                                                  typename GetTypeReal<std::complex<double>>::type* const ekb,
                                                  std::complex<double>* const wfc,
                                                  const MPI_Comm& comm,
                                                  const int diag_subspace,
                                                  const int block_size);

template void diago_hs_para<float>(float* h,
                                   float* s,
                                   const int lda,
                                   const int nband,
                                   typename GetTypeReal<float>::type* const ekb,
                                   float* const wfc,
                                   const MPI_Comm& comm,
                                   const int diag_subspace,
                                   const int block_size);
                                   
template void diago_hs_para<std::complex<float>>(std::complex<float>* h,
                                                 std::complex<float>* s,
                                                 const int lda,
                                                 const int nband,
                                                 typename GetTypeReal<std::complex<float>>::type* const ekb,
                                                 std::complex<float>* const wfc,
                                                 const MPI_Comm& comm,
                                                 const int diag_subspace,
                                                 const int block_size);

}

class TestHSolverPW : public ::testing::Test {
  public:
    // HSolverPW declares this fixture a friend, but a TEST_F body lives in a
    // class derived from it and friendship is not inherited, so the call into
    // the protected hamiltSolvePsiK() is routed through here.
    template <typename T, typename Device>
    static bool hamiltSolvePsiK(hsolver::HSolverPW<T, Device>& hs,
                                const hsolver::HSOperator<T, Device>& op,
                                psi::Psi<T, Device>& ps,
                                std::vector<typename GetTypeReal<T>::type>& pre,
                                typename GetTypeReal<T>::type* eig,
                                const int ntry)
    {
        return hs.hamiltSolvePsiK(op, ps, pre, eig, ntry).converged;
    }

    template <typename T, typename Device>
    static void set_thresholds(hsolver::HSolverPW<T, Device>& hs, const std::vector<double>& thresholds)
    {
        hs.ethr_band = thresholds;
    }

    static void report_failed_kpoints(hsolver::HSolverPW<std::complex<double>>& hs, std::ostream& log)
    {
        using Status = hsolver::HSolverPW<std::complex<double>>::DiagonalizationStatus;
        std::vector<Status> statuses(3);
        statuses[0].converged = false;
        statuses[0].iteration_limit_reached = true;
        statuses[0].active_band_count = 2;
        statuses[0].max_residual = 0.125;
        statuses[2].converged = false;
        statuses[2].active_band_count = 1;
        statuses[2].max_residual = 0.25;
        hs.report_ppcg_convergence(statuses, log);
    }

    ModulePW::PW_Basis_K pwbk;
    hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU> hs_f
        = hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU>(
            &pwbk,
            "scf",
            "pw",
            "cg",
            false,
            1,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::SCF_ITER,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_NMAX,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_THR,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::need_subspace,
            0,
            false,
            4,
            16,
            0,
            0);
    hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hs_d
        = hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU>(
            &pwbk,
            "scf",
            "pw",
            "cg",
            false,
            1,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::SCF_ITER,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_NMAX,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::PW_DIAG_THR,
            hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>::need_subspace,
            0,
            false,
            4,
            16,
            0,
            0);

    IdentityHSOperator<std::complex<double>> hamilt_test_d;
    IdentityHSOperator<std::complex<float>> hamilt_test_f;

    psi::Psi<std::complex<double>> psi_test_cd;
    psi::Psi<std::complex<float>> psi_test_cf;

    elecstate::ElecState elecstate_test;

    std::string method_test = "cg";

    std::vector<float> ekb_f;

    std::ofstream temp_ofs;
};

// TEST_F(TestHSolverPW, solve) {
//     // initial memory and data
//     elecstate_test.ekb.create(1, 2);
//     elecstate_test.pot = new elecstate::Potential;
//     this->ekb_f.resize(2);
//     psi_test_cf.resize(1, 2, 3);
//     psi_test_cd.resize(1, 2, 3);
//     const double nelec = 1.0;

//     // check solve()
//     EXPECT_EQ(this->hs_f.initialed_psi, false);
//     EXPECT_EQ(this->hs_d.initialed_psi, false);

//     this->hs_f.solve(&hamilt_test_f,
//                      psi_test_cf,
//                      &elecstate_test,
//                      elecstate_test.ekb.c,

//                      0,
//                      1,

//                      true);
//     // EXPECT_EQ(this->hs_f.initialed_psi, true);
//     for (int i = 0; i < psi_test_cf.size(); i++) {
//         EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 3);
//     }
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);
//     EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<float>>::avg_iter,
//                      0.0);

//     this->hs_d.solve(&hamilt_test_d,
//                      psi_test_cd,
//                      &elecstate_test,
//                      elecstate_test.ekb.c,

//                      0,
//                      1,

//                      true);
  
//     // EXPECT_EQ(this->hs_d.initialed_psi, true);
//     EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<double>>::avg_iter,
//                      0.0);
//     for (int i = 0; i < psi_test_cd.size(); i++) {
//         EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i + 3);
//     }
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 4.0);
//     EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 7.0);

//     // // check hamiltSolvePsiK()
//     // this->hs_f.hamiltSolvePsiK(&hamilt_test_f, psi_test_cf, this->hs_f.precondition, ekb_f.data());
//     // this->hs_d.hamiltSolvePsiK(&hamilt_test_d,
//     //                            psi_test_cd,
//     //                            this->hs_f.precondition,
//     //                            elecstate_test.ekb.c);
//     // for (int i = 0; i < psi_test_cf.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
//     // }
//     // for (int i = 0; i < psi_test_cd.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
//     // }
//     // EXPECT_DOUBLE_EQ(ekb_f[0], 5.0);
//     // EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[0], 5.0);
//     // EXPECT_DOUBLE_EQ(ekb_f[1], 8.0);
//     // EXPECT_DOUBLE_EQ(elecstate_test.ekb.c[1], 8.0);

//     // // check endDiagH()
//     // this->hs_f.initialed_psi = true;
//     // this->hs_d.initialed_psi = true;
//     // this->hs_f.endDiagh();
//     // this->hs_d.endDiagh();
//     // // will change state of initialed_psi in endDiagh
//     // EXPECT_EQ(this->hs_f.initialed_psi, true);
//     // EXPECT_EQ(this->hs_d.initialed_psi, true);

//     // // check updatePsiK()
//     // // skip initializing Psi, Psi will not change
//     // this->hs_f.updatePsiK(&hamilt_test_f, psi_test_cf, 0);
//     // this->hs_d.updatePsiK(&hamilt_test_d, psi_test_cd, 0);
//     // for (int i = 0; i < psi_test_cf.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cf.get_pointer()[i].real(), i + 4);
//     // }
//     // for (int i = 0; i < psi_test_cd.size(); i++) {
//     //     EXPECT_DOUBLE_EQ(psi_test_cd.get_pointer()[i].real(), i + 4);
//     // }
//     // // check update_precondition()
//     // this->hs_f.update_precondition(this->hs_f.precondition,
//     //                                0,
//     //                                psi_test_cf.get_nbasis());
//     // this->hs_d.update_precondition(this->hs_d.precondition,
//     //                                0,
//     //                                psi_test_cd.get_nbasis());
//     // EXPECT_NEAR(this->hs_f.precondition[0], 2.414213657, 1e-8);
//     // EXPECT_NEAR(this->hs_f.precondition[1], 3.618033886, 1e-8);
//     // EXPECT_NEAR(this->hs_f.precondition[2], 6.236067772, 1e-8);
//     // EXPECT_NEAR(this->hs_d.precondition[0], 2.414213562, 1e-8);
//     // EXPECT_NEAR(this->hs_d.precondition[1], 3.618033989, 1e-8);
//     // EXPECT_NEAR(this->hs_d.precondition[2], 6.236067977, 1e-8);

//     // // check diago_ethr
//     // init_chg = "atomic";
//     // diag_thr = 1e-7;
//     // calculation = "scf";
//     // float test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 1, 1.0);
//     // EXPECT_NEAR(hs_f.diag_ethr, 0.01, 1.0e-7);
//     // EXPECT_NEAR(test_diagethr, 0.01, 1.0e-7);
//     // calculation = "md";
//     // init_chg = "file";
//     // test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 1, 1.0);
//     // EXPECT_NEAR(test_diagethr, 1e-5, 1.0e-7);
//     // test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 2, 1.0);
//     // EXPECT_NEAR(test_diagethr, 0.01, 1.0e-7);
//     // test_diagethr = hs_f.set_diagethr(hs_f.diag_ethr, 0, 3, 1.0e-3);
//     // EXPECT_NEAR(test_diagethr, 0.0001, 1.0e-7);

//     // init_chg = "atomic";
//     // diag_thr = 1e-7;
//     // calculation = "scf";
//     // double test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 1, 1.0);
//     // EXPECT_EQ(hs_d.diag_ethr, 0.01);
//     // EXPECT_EQ(test_diagethr_d, 0.01);
//     // calculation = "md";
//     // init_chg = "file";
//     // test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 1, 1.0);
//     // EXPECT_EQ(test_diagethr_d, 1e-5);
//     // test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 2, 1.0);
//     // EXPECT_EQ(test_diagethr_d, 0.01);
//     // test_diagethr_d = hs_d.set_diagethr(hs_d.diag_ethr, 0, 3, 1.0e-3);
//     // EXPECT_EQ(test_diagethr_d, 0.0001);
// }

TEST_F(TestHSolverPW, SolveLcaoInPW) {
    pwbk.nks = 1;
    // initial memory and data
    elecstate_test.ekb.create(1, 2);
    elecstate_test.wg.create(1,2);
    elecstate_test.klist=new K_Vectors;
    elecstate_test.skip_weights=true;
    elecstate_test.pot = new elecstate::Potential;
    // 1 kpt, 2 bands, 3 basis
    psi_test_cf.resize(1, 2, 3);
    psi_test_cd.resize(1, 2, 3);

    psi::Psi<std::complex<double>> transform_test_cd;
    psi::Psi<std::complex<float>> transform_test_cf;
    // transform psi, the old wanf2, has 2 local basis and 3 pw basis.
    // so in principle the hcc has dimension 3*3 to diagonalize
    // 2 lowest eigenvalues will be selected and save to psi
    transform_test_cd.resize(1, 3, 3);
    transform_test_cf.resize(1, 3, 3);

    // 1, 2, 3 / 4, 5, 6 / 7, 8, 9 would be rank deficient, so the diagonal is
    // lifted to keep the three subspace vectors linearly independent
    for (int iband = 0; iband < transform_test_cd.get_nbands(); iband++) {
        for (int ibasis = 0; ibasis < transform_test_cd.get_nbasis();
             ibasis++) {
            const double value = iband * transform_test_cd.get_nbasis() + ibasis + 1 + (iband == ibasis ? 10.0 : 0.0);
            transform_test_cd
                .get_pointer()[iband * transform_test_cd.get_nbasis() + ibasis]
                = std::complex<double>(value, 0.0);
            transform_test_cf
                .get_pointer()[iband * transform_test_cf.get_nbasis() + ibasis]
                = std::complex<float>(value, 0.0);
        }
    }
    // with H = S = 1 every subspace eigenvalue is 1 and the rotated psi must
    // come out orthonormal
    auto check_orthonormal = [](const auto& p, const double tol) {
        const int nb = p.get_nbands();
        const int nbasis = p.get_nbasis();
        for (int i = 0; i < nb; i++) {
            for (int j = 0; j < nb; j++) {
                std::complex<double> dot = 0.0;
                for (int ig = 0; ig < nbasis; ig++) {
                    dot += std::conj(std::complex<double>(p.get_pointer()[i * nbasis + ig]))
                           * std::complex<double>(p.get_pointer()[j * nbasis + ig]);
                }
                EXPECT_NEAR(dot.real(), i == j ? 1.0 : 0.0, tol);
                EXPECT_NEAR(dot.imag(), 0.0, tol);
            }
        }
    };
    // check solve()
    elecstate_test.ekb.c[0] = 1.0;
    elecstate_test.ekb.c[1] = 2.0;

    hsolver::HSolverLIP<std::complex<float>> hs_f_lip
        = hsolver::HSolverLIP<std::complex<float>>(&pwbk, false, "pw", "scf", elecstate_test.ekb.nc);
    hsolver::HSolverLIP<std::complex<double>> hs_d_lip
        = hsolver::HSolverLIP<std::complex<double>>(&pwbk, false, "pw", "scf", elecstate_test.ekb.nc);
#ifdef __MPI
    const hsolver::diag_comm_info diag_comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info diag_comm(0, 1);
#endif
    std::ostringstream log;
    hs_f_lip.solve(hamilt_test_f,
                   psi_test_cf,
                   &elecstate_test,
                   transform_test_cf,
                   diag_comm,
                   log,
                   true);
    EXPECT_NE(log.str().find("Average iterative diagonalization steps"), std::string::npos);
    EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<float>>::avg_iter, 0.0);
    check_orthonormal(psi_test_cf, 1e-5);
    EXPECT_NEAR(elecstate_test.ekb.c[0], 1.0, 1e-5);
    EXPECT_NEAR(elecstate_test.ekb.c[1], 1.0, 1e-5);

    elecstate_test.ekb.c[0] = 1.0;
    elecstate_test.ekb.c[1] = 2.0;
    hs_d_lip.solve(hamilt_test_d,
                   psi_test_cd,
                   &elecstate_test,
                   transform_test_cd,
                   diag_comm,
                   log,
                   true);
    EXPECT_DOUBLE_EQ(hsolver::DiagoIterAssist<std::complex<double>>::avg_iter, 0.0);
    check_orthonormal(psi_test_cd, 1e-10);
    EXPECT_NEAR(elecstate_test.ekb.c[0], 1.0, 1e-10);
    EXPECT_NEAR(elecstate_test.ekb.c[1], 1.0, 1e-10);
}

// Test that the program exits with an error when npwx < nbands,
// which would cause rank deficiency and psi_norm <= 0 during diagonalization.
TEST_F(TestHSolverPW, NpwxLessThanNbandsDeath)
{
    // Create psi with 5 bands but only 3 basis functions -> npwx=3 < nbands=5
    psi_test_cd.resize(1, 5, 3);
    std::vector<double> precond(3, 0.0);
    std::vector<double> eigenvalues(5, 0.0);
    // Expect death from WARNING_QUIT due to npwx < nbands
    EXPECT_EXIT(
        hamiltSolvePsiK(hs_d, hamilt_test_d, psi_test_cd, precond, eigenvalues.data(), 1),
        ::testing::ExitedWithCode(1),
        ".*"
    );
}

TEST_F(TestHSolverPW, PpcgUsesIndependentConstructorConfiguration)
{
    using T = std::complex<double>;
    const int dim = 18;
    const int nband = 3;
    const int maxiter = 4;
    const double energy_threshold = 1.0e-8;
    const double residual_threshold = 1.0e-4;
    const std::vector<double> energy_thresholds = {1.0e-8, 4.0e-8, 9.0e-8};
    const std::vector<double> residual_thresholds = {1.0e-4, 2.0e-4, 3.0e-4};
    std::vector<double> precondition(dim, 1.0);
    std::vector<T> initial(dim * nband);
    for (size_t i = 0; i < initial.size(); ++i)
    {
        const double index = double(i + 1);
        const double phase = index * 1.7;
        const double real = std::sin(index);
        const double imaginary = std::cos(phase);
        initial[i] = T(real, imaginary);
    }

    // Neither object calls setup_diago_params_pw. Dispatch must honor its own
    // block size and RR interval, including after a different object was used.
    hsolver::HSolverPW<T> first(&pwbk, "scf", "pw", "ppcg", false, 1, 1,
                               maxiter, energy_threshold, false, nband, false, 1, 3, 0, 0);
    hsolver::HSolverPW<T> second(&pwbk, "scf", "pw", "ppcg", false, 1, 1,
                                maxiter, energy_threshold, false, nband, false, 2, 1, 0, 0);
    set_thresholds(first, energy_thresholds);
    set_thresholds(second, energy_thresholds);
    hsolver::HSolverPW<T>* solvers[] = {&first, &second, &first};
    const int block_sizes[] = {1, 2, 1};
    const int rr_steps[] = {3, 1, 3};
    for (int run = 0; run < 3; ++run)
    {
        psi::Psi<T> wavefunctions(1, nband, dim, dim, true);
        wavefunctions.set_all_psi(initial.data(), initial.size());
        std::vector<T> reference = initial;
        std::vector<double> actual_eigenvalues(nband);
        std::vector<double> reference_eigenvalues(nband);
        DiagonalPPCGOperator actual_operator;
        DiagonalPPCGOperator reference_operator;
        hsolver::DiagoPPCG<T, base_device::DEVICE_CPU> direct(
            residual_threshold, maxiter, block_sizes[run], rr_steps[run]);
        const auto hpsi = [&](T* x, T* hx, const int ld, const int nvec) {
            reference_operator.hpsi(x, hx, ld, nvec);
        };
        const auto spsi = [&](T* x, T* sx, const int ld, const int nvec) {
            reference_operator.spsi(x, sx, ld, nvec);
        };
        direct.diag(hpsi, spsi, dim, nband, dim, reference.data(),
                    reference_eigenvalues.data(), residual_thresholds, precondition.data());
        const bool converged = hamiltSolvePsiK(*solvers[run], actual_operator, wavefunctions,
                                               precondition, actual_eigenvalues.data(), 1);
        EXPECT_EQ(converged, direct.converged());
        EXPECT_EQ(actual_operator.block_sizes, reference_operator.block_sizes);
        for (int ib = 0; ib < nband; ++ib)
        {
            EXPECT_NEAR(actual_eigenvalues[ib], reference_eigenvalues[ib], 1.0e-12);
        }
    }
}

TEST_F(TestHSolverPW, PpcgFailureUsesApplicationDiagnostic)
{
    using T = std::complex<double>;
    hsolver::HSolverPW<T> solver(&pwbk, "scf", "pw", "ppcg", false, 1, 1,
                                4, 1.0e-8, false, 1, false, 1, 3, 0, 0);
    const std::vector<double> thresholds(1, 1.0e-8);
    set_thresholds(solver, thresholds);
    psi::Psi<T> wavefunctions(1, 1, 3, 3, true);
    std::vector<double> precondition(3, 1.0);
    precondition[0] = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> eigenvalues(1);
    EXPECT_EXIT({
        std::cout.rdbuf(std::cerr.rdbuf());
        hamiltSolvePsiK(solver, hamilt_test_d, wavefunctions,
                        precondition, eigenvalues.data(), 1);
    }, ::testing::ExitedWithCode(1), "PPCG diagonalization failed:.*preconditioner contains non-finite");
}

TEST_F(TestHSolverPW, PpcgConvergenceWarningSummarizesKpoints)
{
    int contributing_pools = 1;
    Parallel_Reduce::reduce_all(contributing_pools);
    const int kpoints = 2 * contributing_pools;
    const int bands = 3 * contributing_pools;
    const std::string expected_kpoints = std::to_string(kpoints) + " k point(s)";
    const std::string expected_limits = std::to_string(contributing_pools) + " reached the iteration limit";
    const std::string expected_bands = "active bands across those k points = " + std::to_string(bands);
    std::ostringstream log;
    report_failed_kpoints(hs_d, log);
    const std::string message = log.str();
    EXPECT_NE(message.find(expected_kpoints), std::string::npos);
    EXPECT_NE(message.find(expected_limits), std::string::npos);
    EXPECT_NE(message.find(expected_bands), std::string::npos);
    EXPECT_NE(message.find("2.500000e-01"), std::string::npos);
    EXPECT_EQ(std::count(message.begin(), message.end(), '\n'), 1);
}

int main(int argc, char** argv)
{
#ifdef __MPI
    int nproc = 1;
    int myrank = 0;
    setupmpi(argc, argv, nproc, myrank);
    int nproc_in_pool = 1;
    int kpar = nproc;
    int mypool = 0;
    int rank_in_pool = 0;
    divide_pools(nproc, myrank, nproc_in_pool, kpar, mypool, rank_in_pool);
#endif
    testing::InitGoogleTest(&argc, argv);
    // Re-execute death tests so they initialize their own MPI runtime instead
    // of using MPI state inherited from a forked parent process.
    testing::FLAGS_gtest_death_test_style = "threadsafe";
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif
    return result;
}
