#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <streambuf>
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/magnetism.h"
#include "source_cell/pseudo.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "dfpt_stru_fixture.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_global.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"
#include "source_pw/module_pwdft/dftu_base.h"

// ctor/dtor stubs for the cell/spepot link closures live in the shared
// dfpt_test_mocks.cpp compiled into every DFPT test binary.

/************************************************
 *  unit test of DFPT_PW_Data (Phase 4 wiring; B4
 *  absorbed the retired DFPT_IrrepData adapter)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_PW_Data::get_nq() / get_nirr() / get_irrep_modes()
 *     - delegates to the QList irrep data
 *   - DFPT_PW_Data::set/get_dpsi, drho_r, drho_g, dv_r
 *     - per-q storage round trip and out-of-range safety
 *   - DFPT_PW_Data::set/get_converged, add/get_residuals,
 *     set/get_current_iter
 *     - per-(q, irrep) SCF ledger (B4: sunk from DFPT_IrrepData)
 */

class DFPT_PW_DataTest : public DFPTStruTestFixture
{
  protected:
    ModuleCell::QList qlist;
    ModuleDFPT::DFPT_PW_Data data;
    std::ifstream ifs;
    std::ofstream ofs;
    std::ofstream ofs_running;
    std::string output;

    // build a reduced 2x2x2 q-mesh (4 irreducible q-points for O_h)
    void init_qlist()
    {
        construct_ucell(stru_lib[0]);
        ofs_running.open("tmp_dfpt_qlist");
        ModuleSymmetry::Symmetry symm;
        const int cal_symm_repr[2] = {0, 6};
        symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);
        qlist.generate_mesh(ucell, symm, {2, 2, 2}, true);
        data.init(&qlist, 1, 2, 3, 0, 1, 1, nullptr);
    }

    void clear_qlist()
    {
        data.clean();
        ofs_running.close();
        ClearUcell();
        remove("tmp_dfpt_qlist");
    }
};

TEST_F(DFPT_PW_DataTest, DelegatesToQList)
{
    init_qlist();

    EXPECT_EQ(data.get_nq(), qlist.get_nq());
    EXPECT_EQ(data.get_nq(), 4);
    for (int q_idx = 0; q_idx < data.get_nq(); ++q_idx)
    {
        EXPECT_EQ(data.get_nirr(q_idx), 1);
        EXPECT_TRUE(data.get_irrep_modes(q_idx, 0).empty());
    }

    // the first irreducible q-point must be Gamma
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).y, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).z, 0.0);

    clear_qlist();
}

TEST_F(DFPT_PW_DataTest, AccessorsAreBoundSafe)
{
    init_qlist();

    const int nq = data.get_nq();

    // out-of-range access must be safe and return empty containers
    EXPECT_TRUE(data.get_irrep_modes(-1, 0).empty());
    EXPECT_TRUE(data.get_irrep_modes(nq, 0).empty());
    EXPECT_TRUE(data.get_dpsi(-1, 0, 0).empty());
    EXPECT_TRUE(data.get_drho_r(0, 5).empty());
    EXPECT_TRUE(data.get_drho_g(0, 5).empty());
    EXPECT_TRUE(data.get_dv_r(0, 5).empty());

    clear_qlist();
}

TEST_F(DFPT_PW_DataTest, SetterRoundTrip)
{
    init_qlist();

    // dpsi / drho / dv storage went live with C1 (dv) and C3 (drho):
    // reads must return what was written through the per-q slots
    std::vector<std::complex<double>> psi(3, std::complex<double>(1.0, 2.0));
    data.set_dpsi(0, 0, 0, psi);
    std::vector<double> rho(2, 3.0);
    data.set_drho_r(0, 0, rho);
    data.set_drho_g(0, 0, std::vector<std::complex<double>>(2, std::complex<double>(1.0, 0.0)));
    data.set_dv_r(0, 0, rho);

    EXPECT_FALSE(data.get_dpsi(0, 0, 0).empty());
    EXPECT_FALSE(data.get_drho_r(0, 0).empty());
    EXPECT_FALSE(data.get_drho_g(0, 0).empty());
    EXPECT_FALSE(data.get_dv_r(0, 0).empty());

    clear_qlist();
}

TEST_F(DFPT_PW_DataTest, PerIrrepScfLedger)
{
    init_qlist();

    // the ledger must be independent per (q_idx, irrep) — the shape DFPT_PW
    // ::run drives and the stage-A irrep implementation will fill (B4:
    // absorbed from the retired DFPT_IrrepData adapter)
    data.set_converged(0, 0, false);
    data.set_converged(1, 0, true);
    EXPECT_FALSE(data.get_converged(0, 0));
    EXPECT_TRUE(data.get_converged(1, 0));

    data.add_residual(0, 0, 1e-3);
    data.add_residual(0, 0, 2e-4);
    data.add_residual(1, 0, 9e-5);
    EXPECT_EQ(data.get_residuals(0, 0).size(), 2);
    EXPECT_EQ(data.get_residuals(1, 0).size(), 1);
    EXPECT_NEAR(data.get_residuals(0, 0)[1], 2e-4, 1e-12);

    data.set_current_iter(0, 0, 3);
    EXPECT_EQ(data.get_current_iter(0, 0), 3);
    EXPECT_EQ(data.get_current_iter(1, 0), 0); // untouched key defaults to 0
    EXPECT_FALSE(data.get_converged(2, 0));    // untouched key defaults to false
    EXPECT_TRUE(data.get_residuals(2, 0).empty());

    // clean() drops the whole ledger
    data.clean();
    data.init(&qlist, 1, 2, 3, 0, 1, 1, nullptr);
    EXPECT_TRUE(data.get_residuals(0, 0).empty());
    EXPECT_EQ(data.get_current_iter(0, 0), 0);
    EXPECT_FALSE(data.get_converged(1, 0));

    clear_qlist();
}

TEST_F(DFPT_PW_DataTest, DftuReservationWithNullProvider)
{
    init_qlist();

    // no Plus_U wired -> with_u()/u_active() must both be false and docc
    // reads must return empty, so the no-DFT+U path is untouched (U0)
    EXPECT_FALSE(data.with_u());
    EXPECT_FALSE(data.u_active());
    EXPECT_EQ(data.get_dftu(), nullptr);
    EXPECT_TRUE(data.get_docc(0).empty());

    // docc storage is independent of the provider: roundtrip works even
    // with a null provider, out-of-range reads stay safe (U0)
    std::vector<std::complex<double>> occ(4, std::complex<double>(0.5, 0.0));
    data.set_docc(0, occ);
    data.set_docc(2, occ);
    ASSERT_EQ(data.get_docc(0).size(), 4);
    EXPECT_DOUBLE_EQ(data.get_docc(0)[1].real(), 0.5);
    EXPECT_DOUBLE_EQ(data.get_docc(0)[1].imag(), 0.0);
    EXPECT_TRUE(data.get_docc(1).empty());
    EXPECT_TRUE(data.get_docc(-1).empty());
    EXPECT_TRUE(data.get_docc(7).empty());

    clear_qlist();
}

TEST_F(DFPT_PW_DataTest, DftuReservationProviderUsability)
{
    init_qlist();

    // U0 reservation semantics after the Plus_U_Base migration: with_u()
    // reports the wiring, u_active() follows the occupation-matrix state of
    // the provider. A default-constructed provider (occupation matrices not
    // initialized, e.g. no DFT+U ran in the ground state) stays inactive;
    // once the provider marks its occupation matrices initialized (a real
    // ground-state DFT+U run), the reservation turns active. DFPT_PW::init
    // rejects a wired provider outright; this pins the data-level contract.
    Plus_U_Base dftu;
    data.init(&qlist, 1, 2, 3, 0, 1, 1, &dftu);
    EXPECT_TRUE(data.with_u());
    EXPECT_FALSE(data.u_active());
    ASSERT_NE(data.get_dftu(), nullptr);

    dftu.mark_occ_mat_initialized();
    EXPECT_TRUE(data.u_active());

    data.clean();
    // clean() resets the response storage but keeps the provider wiring
    EXPECT_TRUE(data.with_u());

    clear_qlist();
}
