#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <iostream>
#include <streambuf>
#define private public
#include "source_cell/atom_pseudo.h"
#include "source_cell/atom_spec.h"
#include "source_cell/pseudo.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_cell/magnetism.h"
#undef private
#include "source_base/mathzone.h"
#include "source_base/parallel_global.h"
#include "source_base/global_variable.h"
#include "source_pw/module_dfpt/dfpt_irrep_data.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"

pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
Atom::Atom()
{
}
Atom::~Atom()
{
}
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
SepPot::SepPot() {}
SepPot::~SepPot() {}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

/************************************************
 *  unit test of DFPT_IrrepData (Phase 4 wiring)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_IrrepData::get_nq() / get_nirr() / get_irrep_modes()
 *     - delegates to the underlying QList irrep data
 *   - DFPT_IrrepData::set/get_dpsi, drho_r, drho_g, dv_r
 *     - irrep-indexed accessors forward to the per-q DFPT_PW_Data
 *   - DFPT_IrrepData::set/get_converged, add/get_residuals,
 *     set/get_current_iter
 *     - per-(q, irrep) SCF bookkeeping
 */

// abbreviated from module_symmetry/test/symm_test.cpp and klist_test.cpp
struct atomtype_
{
    std::string atomname;
    std::vector<std::vector<double>> coordinate;
};

struct stru_
{
    int ibrav;
    std::string point_group;    // Schoenflies symbol
    std::string point_group_hm; // Hermann-Mauguin notation.
    std::string space_group;
    std::vector<double> cell;
    std::vector<atomtype_> all_type;
};

std::vector<stru_> stru_lib{stru_{1,
                                  "O_h",
                                  "m-3m",
                                  "Pm-3m",
                                  std::vector<double>{1., 0., 0., 0., 1., 0., 0., 0., 1.},
                                  std::vector<atomtype_>{atomtype_{"C",
                                                                   std::vector<std::vector<double>>{
                                                                       {0., 0., 0.},
                                                                   }}}}};

class DFPT_IrrepDataTest : public testing::Test
{
  protected:
    ModuleCell::QList qlist;
    ModuleDFPT::DFPT_PW_Data data;
    std::ifstream ifs;
    std::ofstream ofs;
    std::ofstream ofs_running;
    std::string output;

    UnitCell ucell;
    void construct_ucell(stru_& stru)
    {
        std::vector<atomtype_> coord = stru.all_type;
        ucell.a1 = ModuleBase::Vector3<double>(stru.cell[0], stru.cell[1], stru.cell[2]);
        ucell.a2 = ModuleBase::Vector3<double>(stru.cell[3], stru.cell[4], stru.cell[5]);
        ucell.a3 = ModuleBase::Vector3<double>(stru.cell[6], stru.cell[7], stru.cell[8]);
        ucell.ntype = stru.all_type.size();
        ucell.atoms = new Atom[ucell.ntype];
        ucell.nat = 0;
        ucell.latvec.e11 = ucell.a1.x;
        ucell.latvec.e12 = ucell.a1.y;
        ucell.latvec.e13 = ucell.a1.z;
        ucell.latvec.e21 = ucell.a2.x;
        ucell.latvec.e22 = ucell.a2.y;
        ucell.latvec.e23 = ucell.a2.z;
        ucell.latvec.e31 = ucell.a3.x;
        ucell.latvec.e32 = ucell.a3.y;
        ucell.latvec.e33 = ucell.a3.z;
        ucell.GT = ucell.latvec.Inverse();
        ucell.G = ucell.GT.Transpose();
        ucell.lat0 = 1.8897261254578281;
        for (int i = 0; i < coord.size(); i++)
        {
            ucell.atoms[i].label = coord[i].atomname;
            ucell.atoms[i].na = coord[i].coordinate.size();
            ucell.atoms[i].tau.resize(ucell.atoms[i].na);
            ucell.atoms[i].taud.resize(ucell.atoms[i].na);
            for (int j = 0; j < ucell.atoms[i].na; j++)
            {
                std::vector<double> this_atom = coord[i].coordinate[j];
                ucell.atoms[i].tau[j] = ModuleBase::Vector3<double>(this_atom[0], this_atom[1], this_atom[2]);
                ModuleBase::Mathzone::Cartesian_to_Direct(ucell.atoms[i].tau[j].x,
                                                          ucell.atoms[i].tau[j].y,
                                                          ucell.atoms[i].tau[j].z,
                                                          ucell.a1.x,
                                                          ucell.a1.y,
                                                          ucell.a1.z,
                                                          ucell.a2.x,
                                                          ucell.a2.y,
                                                          ucell.a2.z,
                                                          ucell.a3.x,
                                                          ucell.a3.y,
                                                          ucell.a3.z,
                                                          ucell.atoms[i].taud[j].x,
                                                          ucell.atoms[i].taud[j].y,
                                                          ucell.atoms[i].taud[j].z);
            }
            ucell.nat += ucell.atoms[i].na;
        }
    }

    void ClearUcell()
    {
        delete[] ucell.atoms;
    }

    // build a reduced 2x2x2 q-mesh (4 irreducible q-points for O_h)
    void init_qlist()
    {
        construct_ucell(stru_lib[0]);
        GlobalV::ofs_running.open("tmp_dfpt_qlist");
        ModuleSymmetry::Symmetry symm;
        const int cal_symm_repr[2] = {0, 6};
        symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, GlobalV::ofs_running, 1e-6, 1, "scf", cal_symm_repr);
        qlist.generate_mesh(ucell, symm, {2, 2, 2}, true);
        data.init(&qlist, 1, 2, 3, 0, 1, 1, nullptr);
    }

    void clear_qlist()
    {
        data.clean();
        GlobalV::ofs_running.close();
        ClearUcell();
        remove("tmp_dfpt_qlist");
    }
};

TEST_F(DFPT_IrrepDataTest, DelegatesToQList)
{
    init_qlist();

    ModuleDFPT::DFPT_IrrepData irrep_data(data);

    EXPECT_EQ(irrep_data.get_nq(), qlist.get_nq());
    EXPECT_EQ(irrep_data.get_nq(), 4);
    for (int q_idx = 0; q_idx < irrep_data.get_nq(); ++q_idx)
    {
        EXPECT_EQ(irrep_data.get_nirr(q_idx), 1);
        EXPECT_TRUE(irrep_data.get_irrep_modes(q_idx, 0).empty());
    }

    // the first irreducible q-point must be Gamma
    EXPECT_DOUBLE_EQ(qlist.get_q(0).x, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).y, 0.0);
    EXPECT_DOUBLE_EQ(qlist.get_q(0).z, 0.0);

    clear_qlist();
}

TEST_F(DFPT_IrrepDataTest, IrrepIndexedAccessorsAreBoundSafe)
{
    init_qlist();

    ModuleDFPT::DFPT_IrrepData irrep_data(data);
    const int nq = irrep_data.get_nq();
    const int nirr = irrep_data.get_nirr(0);

    // out-of-range access must be safe and return empty containers
    EXPECT_TRUE(irrep_data.get_irrep_modes(-1, 0).empty());
    EXPECT_TRUE(irrep_data.get_irrep_modes(nq, 0).empty());
    EXPECT_TRUE(irrep_data.get_dpsi(-1, 0, 0, 0).empty());
    EXPECT_TRUE(irrep_data.get_drho_r(0, 5, 0).empty());
    EXPECT_TRUE(irrep_data.get_drho_g(0, 5, 0).empty());
    EXPECT_TRUE(irrep_data.get_dv_r(0, 5, 0).empty());

    (void)nirr;

    clear_qlist();
}

TEST_F(DFPT_IrrepDataTest, SetterRoundTripViaWrapper)
{
    init_qlist();

    ModuleDFPT::DFPT_IrrepData irrep_data(data);

    // dpsi / drho / dv are still design-phase no-op storage stubs: the
    // wrapper must forward the irrep-indexed calls to the per-q storage
    // slot without crashing, and reads must reflect the stub (empty)
    std::vector<std::complex<double>> psi(3, std::complex<double>(1.0, 2.0));
    irrep_data.set_dpsi(0, 0, 0, 0, psi);
    std::vector<double> rho(2, 3.0);
    irrep_data.set_drho_r(0, 0, 0, rho);
    irrep_data.set_drho_g(0, 0, 0, std::vector<std::complex<double>>(2, std::complex<double>(1.0, 0.0)));
    irrep_data.set_dv_r(0, 0, 0, rho);

    // the wrapper reads the same slot the setter wrote through; slots that
    // are now backed by real storage return non-empty, while design-phase
    // stubs still return empty.
    EXPECT_FALSE(irrep_data.get_dpsi(0, 0, 0, 0).empty());
    EXPECT_TRUE(irrep_data.get_drho_r(0, 0, 0).empty());
    EXPECT_TRUE(irrep_data.get_drho_g(0, 0, 0).empty());
    EXPECT_FALSE(irrep_data.get_dv_r(0, 0, 0).empty());

    clear_qlist();
}

TEST_F(DFPT_IrrepDataTest, PerIrrepScfBookkeeping)
{
    init_qlist();

    ModuleDFPT::DFPT_IrrepData irrep_data(data);
    const int nirr = irrep_data.get_nirr(0);

    // bookkeeping must be independent per (q_idx, irrep)
    irrep_data.set_converged(0, 0, false);
    irrep_data.set_converged(1, 0, true);
    EXPECT_FALSE(irrep_data.get_converged(0, 0));
    EXPECT_TRUE(irrep_data.get_converged(1, 0));

    irrep_data.add_residual(0, 0, 1e-3);
    irrep_data.add_residual(0, 0, 2e-4);
    irrep_data.add_residual(1, 0, 9e-5);
    EXPECT_EQ(irrep_data.get_residuals(0, 0).size(), 2);
    EXPECT_EQ(irrep_data.get_residuals(1, 0).size(), 1);
    EXPECT_NEAR(irrep_data.get_residuals(0, 0)[1], 2e-4, 1e-12);

    irrep_data.set_current_iter(0, 0, 3);
    EXPECT_EQ(irrep_data.get_current_iter(0, 0), 3);
    EXPECT_EQ(irrep_data.get_current_iter(1, 0), 0); // untouched key defaults to 0

    for (int irrep = 0; irrep < nirr; ++irrep)
    {
        EXPECT_FALSE(irrep_data.get_converged(0, irrep));
    }

    clear_qlist();
}

TEST_F(DFPT_IrrepDataTest, DftuReservationWithNullProvider)
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
