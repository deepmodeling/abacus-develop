#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>

#include "hsolver_pw_sup.h"
#include "hsolver_supplementary_mock.h"
#include "source_hsolver/diag_comm_info.h"
#include "source_hsolver/diago_params.h"
#include "source_hsolver/hs_operator.h"
#include "source_hsolver/hsolver_lcaopw.h"
#include "source_hsolver/hsolver_pw.h"
#ifdef __MPI
#include "source_basis/module_pw/test/test_tool.h"
#endif

#include <algorithm>
#include <cmath>
#include <complex>

/// H = S = identity: the simplest operator the solvers can be handed
template <typename T> class IdentityHSOperator : public hsolver::HSOperator<T, base_device::DEVICE_CPU>
{
  public:
    void update_k(const int ik) override {}
    void hpsi(
        const T *x,
        T *hx,
        const int ld,
        const int nvec) const override
    {
        std::copy(x, x + static_cast<size_t>(ld) * nvec, hx);
    }
    void spsi(
        const T *x,
        T *sx,
        const int ld,
        const int nvec) const override
    {
        std::copy(x, x + static_cast<size_t>(ld) * nvec, sx);
    }
};

class DiagonalHSOperator : public hsolver::HSOperator<std::complex<double>, base_device::DEVICE_CPU>
{
  public:
    explicit DiagonalHSOperator(const std::vector<double>& diagonal) : diagonal_(diagonal) {}

    void update_k(const int ik) override {}

    void hpsi(const std::complex<double>* x,
              std::complex<double>* hx,
              const int ld,
              const int nvec) const override
    {
        for (int column = 0; column < nvec; ++column)
        {
            for (int row = 0; row < ld; ++row)
            {
                hx[row + column * ld] = diagonal_[row] * x[row + column * ld];
            }
        }
    }

    void spsi(const std::complex<double>* x,
              std::complex<double>* sx,
              const int ld,
              const int nvec) const override
    {
        const size_t count = static_cast<size_t>(ld) * nvec;
        const std::complex<double>* end = x + count;
        std::copy(x, end, sx);
    }

  private:
    std::vector<double> diagonal_;
};

// Mock implementations for the template functions causing linking errors
namespace ModulePW
{
// Mock implementation for recip_to_real
template <typename FPTYPE, typename Device>
void PW_Basis_K::recip_to_real(
    const Device *ctx,
    const std::complex<FPTYPE> *in,
    std::complex<FPTYPE> *out,
    const int ik,
    const bool add,
    const FPTYPE factor) const
{
    // Simple mock implementation that does nothing
    // In a real test, you might want to implement behavior that simulates the actual function
}

// Mock implementation for real_to_recip
template <typename FPTYPE, typename Device>
void PW_Basis_K::real_to_recip(
    const Device *ctx,
    const std::complex<FPTYPE> *in,
    std::complex<FPTYPE> *out,
    const int ik,
    const bool add,
    const FPTYPE factor) const
{
    // Simple mock implementation that does nothing
}

// Explicit template instantiations
template void PW_Basis_K::recip_to_real<float, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU *ctx,
    const std::complex<float> *in,
    std::complex<float> *out,
    const int ik,
    const bool add,
    const float factor) const;

template void PW_Basis_K::recip_to_real<double, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU *ctx,
    const std::complex<double> *in,
    std::complex<double> *out,
    const int ik,
    const bool add,
    const double factor) const;

template void PW_Basis_K::real_to_recip<float, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU *ctx,
    const std::complex<float> *in,
    std::complex<float> *out,
    const int ik,
    const bool add,
    const float factor) const;

template void PW_Basis_K::real_to_recip<double, base_device::DEVICE_CPU>(
    const base_device::DEVICE_CPU *ctx,
    const std::complex<double> *in,
    std::complex<double> *out,
    const int ik,
    const bool add,
    const double factor) const;
} // namespace ModulePW

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
namespace hsolver
{
template <typename T>
void diago_hs_para(
    T *h,
    T *s,
    const int lda,
    const int nband,
    typename GetTypeReal<T>::type *const ekb,
    T *const wfc,
    const MPI_Comm &comm,
    const int diag_subspace,
    const int block_size = 0)
{
}
template void diago_hs_para<double>(
    double *h,
    double *s,
    const int lda,
    const int nband,
    typename GetTypeReal<double>::type *const ekb,
    double *const wfc,
    const MPI_Comm &comm,
    const int diag_subspace,
    const int block_size);

template void diago_hs_para<std::complex<double>>(
    std::complex<double> *h,
    std::complex<double> *s,
    const int lda,
    const int nband,
    typename GetTypeReal<std::complex<double>>::type *const ekb,
    std::complex<double> *const wfc,
    const MPI_Comm &comm,
    const int diag_subspace,
    const int block_size);

template void diago_hs_para<float>(
    float *h,
    float *s,
    const int lda,
    const int nband,
    typename GetTypeReal<float>::type *const ekb,
    float *const wfc,
    const MPI_Comm &comm,
    const int diag_subspace,
    const int block_size);

template void diago_hs_para<std::complex<float>>(
    std::complex<float> *h,
    std::complex<float> *s,
    const int lda,
    const int nband,
    typename GetTypeReal<std::complex<float>>::type *const ekb,
    std::complex<float> *const wfc,
    const MPI_Comm &comm,
    const int diag_subspace,
    const int block_size);

} // namespace hsolver

class TestHSolverPW : public ::testing::Test
{
  public:
    using DiagoAssist = hsolver::DiagoIterAssist<std::complex<double>, base_device::DEVICE_CPU>;

    void SetUp() override
    {
        saved_diag_iter_max = DiagoAssist::PW_DIAG_NMAX;
        saved_diag_thr = DiagoAssist::PW_DIAG_THR;
        saved_scf_iter = DiagoAssist::SCF_ITER;
        saved_need_subspace = DiagoAssist::need_subspace;
        saved_avg_iter = DiagoAssist::avg_iter;
    }

    void TearDown() override
    {
        DiagoAssist::PW_DIAG_NMAX = saved_diag_iter_max;
        DiagoAssist::PW_DIAG_THR = saved_diag_thr;
        DiagoAssist::SCF_ITER = saved_scf_iter;
        DiagoAssist::need_subspace = saved_need_subspace;
        DiagoAssist::avg_iter = saved_avg_iter;
    }

    // HSolverPW declares this fixture a friend, but a TEST_F body lives in a
    // class derived from it and friendship is not inherited, so the call into
    // the protected hamiltSolvePsiK() is routed through here.
    template <typename T, typename Device>
    static void hamiltSolvePsiK(
        hsolver::HSolverPW<T, Device> &hs,
        const hsolver::HSOperator<T, Device> &op,
        psi::Psi<T, Device> &ps,
        std::vector<typename GetTypeReal<T>::type> &pre,
        typename GetTypeReal<T>::type *eig,
        const int ntry,
        std::ostream& log)
    {
        hs.hamiltSolvePsiK(
            op,
            ps,
            pre,
            eig,
            ntry,
            log);
    }

    template <typename T, typename Device>
    static void set_ethr_band(hsolver::HSolverPW<T, Device>& hs, const int count, const double value)
    {
        hs.ethr_band.assign(count, value);
    }

    std::vector<double> run_ppcg_dispatch(const std::string& calculation,
                                         const bool exact_vectors,
                                         const int iteration_budget,
                                         std::ostream& log)
    {
        Input_para inp;
        inp.calculation = calculation;
        inp.ks_solver = "ppcg";
        inp.pw_diag_nmax = iteration_budget;
        hsolver::setup_diago_params_pw<std::complex<double>, base_device::DEVICE_CPU>(0, 1, 1.0e-12, inp);
        hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hs(
            &pwbk,
            calculation,
            "pw",
            "ppcg",
            false,
            1,
            DiagoAssist::SCF_ITER,
            DiagoAssist::PW_DIAG_NMAX,
            DiagoAssist::PW_DIAG_THR,
            DiagoAssist::need_subspace,
            2,
            false,
            1,
            0,
            0);
        EXPECT_EQ(hs.diag_iter_max, iteration_budget);
        const std::vector<double> diagonal = {1.0, 2.0, 4.0};
        DiagonalHSOperator op(diagonal);
        psi::Psi<std::complex<double>> ps;
        ps.resize(1, 2, 3);
        std::complex<double>* vectors = ps.get_pointer();
        std::complex<double>* end = vectors + 6;
        const std::complex<double> zero(0.0, 0.0);
        std::fill(vectors, end, zero);
        vectors[0] = std::complex<double>(1.0, 0.0);
        vectors[4] = std::complex<double>(1.0, 0.0);
        if (!exact_vectors)
        {
            vectors[1] = std::complex<double>(1.0, 0.0);
            vectors[5] = std::complex<double>(1.0, 0.0);
        }
        std::vector<double> preconditioner(3, 1.0);
        std::vector<double> eigenvalues(2, 0.0);
        set_ethr_band(hs, 2, 1.0e-12);
        hamiltSolvePsiK(hs, op, ps, preconditioner, eigenvalues.data(), 1, log);
        return eigenvalues;
    }

    ModulePW::PW_Basis_K pwbk;
    hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU> hs_f =
        hsolver::HSolverPW<std::complex<float>, base_device::DEVICE_CPU>(
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
            0,
            0);
    hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU> hs_d =
        hsolver::HSolverPW<std::complex<double>, base_device::DEVICE_CPU>(
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

  private:
    int saved_diag_iter_max = 0;
    double saved_diag_thr = 0.0;
    int saved_scf_iter = 0;
    bool saved_need_subspace = false;
    double saved_avg_iter = 0.0;
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

TEST_F(TestHSolverPW, SolveLcaoInPW)
{
    pwbk.nks = 1;
    // initial memory and data
    elecstate_test.ekb.create(1, 2);
    elecstate_test.wg.create(1, 2);
    elecstate_test.klist = new K_Vectors;
    elecstate_test.skip_weights = true;
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
    for (int iband = 0; iband < transform_test_cd.get_nbands(); iband++)
    {
        for (int ibasis = 0; ibasis < transform_test_cd.get_nbasis(); ibasis++)
        {
            const double value = iband * transform_test_cd.get_nbasis() + ibasis + 1 + (iband == ibasis ? 10.0 : 0.0);
            transform_test_cd.get_pointer()[iband * transform_test_cd.get_nbasis() + ibasis] =
                std::complex<double>(value, 0.0);
            transform_test_cf.get_pointer()[iband * transform_test_cf.get_nbasis() + ibasis] =
                std::complex<float>(value, 0.0);
        }
    }
    // with H = S = 1 every subspace eigenvalue is 1 and the rotated psi must
    // come out orthonormal
    auto check_orthonormal = [](const auto &p, const double tol)
    {
        const int nb = p.get_nbands();
        const int nbasis = p.get_nbasis();
        for (int i = 0; i < nb; i++)
        {
            for (int j = 0; j < nb; j++)
            {
                std::complex<double> dot = 0.0;
                for (int ig = 0; ig < nbasis; ig++)
                {
                    dot += std::conj(std::complex<double>(p.get_pointer()[i * nbasis + ig])) *
                           std::complex<double>(p.get_pointer()[j * nbasis + ig]);
                }
                EXPECT_NEAR(dot.real(), i == j ? 1.0 : 0.0, tol);
                EXPECT_NEAR(dot.imag(), 0.0, tol);
            }
        }
    };
    // check solve()
    elecstate_test.ekb.c[0] = 1.0;
    elecstate_test.ekb.c[1] = 2.0;

    hsolver::HSolverLIP<std::complex<float>> hs_f_lip = hsolver::HSolverLIP<std::complex<float>>(
        &pwbk,
        false,
        "pw",
        "scf",
        elecstate_test.ekb.nc);
    hsolver::HSolverLIP<std::complex<double>> hs_d_lip = hsolver::HSolverLIP<std::complex<double>>(
        &pwbk,
        false,
        "pw",
        "scf",
        elecstate_test.ekb.nc);
#ifdef __MPI
    const hsolver::diag_comm_info diag_comm(MPI_COMM_SELF, 0, 1);
#else
    const hsolver::diag_comm_info diag_comm(0, 1);
#endif
    std::ostringstream log;
    hs_f_lip.solve(
        hamilt_test_f,
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
    hs_d_lip.solve(
        hamilt_test_d,
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
    std::ostringstream log;
    // Expect death from WARNING_QUIT due to npwx < nbands
    EXPECT_EXIT(
        hamiltSolvePsiK(
            hs_d,
            hamilt_test_d,
            psi_test_cd,
            precond,
            eigenvalues.data(),
            1,
            log),
        ::testing::ExitedWithCode(1), ".*");
}

TEST_F(TestHSolverPW, PPCGScfDispatchReportsExhaustedBudget)
{
    std::ostringstream log;
    const std::vector<double> eigenvalues = run_ppcg_dispatch("scf", false, 0, log);
    EXPECT_NE(log.str().find("PPCG did not converge"), std::string::npos);
    EXPECT_NE(log.str().find("active bands=2"), std::string::npos);
    EXPECT_NE(log.str().find("max residual="), std::string::npos);
    EXPECT_NE(log.str().find("iteration limit reached"), std::string::npos);
    for (const double eigenvalue : eigenvalues)
    {
        EXPECT_TRUE(std::isfinite(eigenvalue));
    }
}

TEST_F(TestHSolverPW, PPCGNscfDispatchRejectsExhaustedBudget)
{
    std::ostringstream log;
    EXPECT_EXIT(
        run_ppcg_dispatch("nscf", false, 0, log),
        ::testing::ExitedWithCode(1),
        "PPCG did not converge.*iteration limit reached");
}

TEST_F(TestHSolverPW, PPCGConvergedDispatchHasNoWarning)
{
    std::ostringstream log;
    const std::vector<double> eigenvalues = run_ppcg_dispatch("scf", true, 0, log);
    EXPECT_EQ(log.str().find("PPCG did not converge"), std::string::npos);
    EXPECT_NEAR(eigenvalues[0], 1.0, 1.0e-12);
    EXPECT_NEAR(eigenvalues[1], 2.0, 1.0e-12);
}

TEST_F(TestHSolverPW, PPCGNscfDispatchUsesConfiguredBudget)
{
    DiagoAssist::PW_DIAG_NMAX = 30;
    std::ostringstream log;
    const std::vector<double> eigenvalues = run_ppcg_dispatch("nscf", false, 77, log);
    EXPECT_EQ(log.str().find("PPCG did not converge"), std::string::npos);
    EXPECT_NEAR(eigenvalues[0], 1.0, 1.0e-12);
    EXPECT_NEAR(eigenvalues[1], 2.0, 1.0e-12);
}

TEST_F(TestHSolverPW, DiagoParamsPreservesOtherNscfSolverBudgets)
{
    Input_para inp;
    inp.calculation = "nscf";
    inp.pw_diag_nmax = 77;
    for (const std::string solver : {"cg", "dav", "dav_subspace", "bpcg"})
    {
        SCOPED_TRACE(solver);
        inp.ks_solver = solver;
        DiagoAssist::PW_DIAG_NMAX = 30;
        hsolver::setup_diago_params_pw<std::complex<double>, base_device::DEVICE_CPU>(0, 1, 1.0e-12, inp);
        EXPECT_EQ(DiagoAssist::PW_DIAG_NMAX, 30);
    }
}

TEST_F(TestHSolverPW, DiagoParamsSetsScfSolverBudgets)
{
    Input_para inp;
    inp.calculation = "scf";
    inp.pw_diag_nmax = 77;
    for (const std::string solver : {"cg", "dav", "dav_subspace", "bpcg", "ppcg"})
    {
        SCOPED_TRACE(solver);
        inp.ks_solver = solver;
        DiagoAssist::PW_DIAG_NMAX = 30;
        hsolver::setup_diago_params_pw<std::complex<double>, base_device::DEVICE_CPU>(0, 1, 1.0e-12, inp);
        EXPECT_EQ(DiagoAssist::PW_DIAG_NMAX, 77);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
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
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    finishmpi();
#endif
    return result;
}
