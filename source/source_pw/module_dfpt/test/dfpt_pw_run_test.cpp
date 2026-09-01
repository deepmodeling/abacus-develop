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
#include "source_estate/module_charge/charge_mixing.h"
#include "source_pw/module_pwdft/dftu_base.h"
#include "source_pw/module_dfpt/dfpt_pw.h"

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

Charge_Mixing::~Charge_Mixing() {}

/************************************************
 *  unit test of DFPT_PW::run() (Phase 4 wiring)
 ***********************************************/

/**
 * - Tested Functions:
 *   - DFPT_PW::init() with a reduced q-mesh
 *   - DFPT_PW::run() - the per-irrep SCF loop skeleton. All heavy
 *     solvers (pert / stern / rho) are still design-phase stubs, so
 *     the loop must converge on the first iteration and still invoke
 *     the phonon assembly / diagonalization for every irreducible q.
 *   - DFPT_PW::get_phonon_freq() - one frequency per phonon mode,
 *     i.e. 3*nat entries for each q.
 */

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

class DFPT_PWRunTest : public testing::Test
{
  protected:
    ModuleDFPT::DFPT_PW dfpt;
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

    void SetUp() override
    {
        construct_ucell(stru_lib[0]);
        ofs_running.open("tmp_dfpt_run");
        ModuleSymmetry::Symmetry symm;
        const int cal_symm_repr[2] = {0, 6};
        symm.analy_sys(ucell.lat, ucell.st, ucell.atoms, ofs_running, 1e-6, 1, "scf", cal_symm_repr);
        ucell.symm = symm;
    }

    void TearDown() override
    {
        ofs_running.close();
        ClearUcell();
        remove("tmp_dfpt_run");
    }
};

TEST_F(DFPT_PWRunTest, RunsPerIrrepLoopForAllQ)
{
    dfpt.set_qmesh(2, 2, 2); // reduced to 4 irreducible q in O_h
    dfpt.set_max_iter(10);
    psi::Psi<std::complex<double>> psi;
    // skeleton mode: no bases wired (design-phase fallback of the irrep loop)
    dfpt.init(ucell, psi, nullptr, nullptr, nullptr, std::vector<double>(),
              ModuleBase::matrix(), ModuleBase::matrix(), nullptr, 1.0, 15.0, nullptr);
    dfpt.run();

    // each of the 4 irreducible q points must expose 3*nat phonon modes
    const int expected_modes = 3 * ucell.nat;
    for (int q_idx = 0; q_idx < 4; ++q_idx)
    {
        EXPECT_EQ(dfpt.get_phonon_freq(q_idx).size(), expected_modes);
    }
}

TEST_F(DFPT_PWRunTest, DielectricAndBornAreExposed)
{
    dfpt.set_qmesh(1, 1, 1); // Gamma-only q mesh
    psi::Psi<std::complex<double>> psi;
    dfpt.init(ucell, psi, nullptr, nullptr, nullptr, std::vector<double>(),
              ModuleBase::matrix(), ModuleBase::matrix(), nullptr, 1.0, 15.0, nullptr);
    dfpt.run();

    // design-phase stubs return default-constructed matrices
    ModuleBase::matrix eps = dfpt.get_dielectric_tensor();
    ModuleBase::matrix born = dfpt.get_born_charges(0);
    EXPECT_EQ(eps.nr, 0);
    EXPECT_EQ(eps.nc, 0);
    EXPECT_EQ(born.nr, 0);
    EXPECT_EQ(born.nc, 0);
}

TEST_F(DFPT_PWRunTest, DftuReservationWithProviderRejectsInit)
{
    // the preceding irrep-loop tests run OpenMP regions, so the default
    // "fast" fork-based death test can deadlock in the multithreaded child
    // (CI: OMP_NUM_THREADS=2); use the fork+exec style instead
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    // DFT+U reservation (U0): the ground state wires a provider when
    // dft_plus_u is enabled upstream, and PW-basis DFT+U actually runs in
    // the ground state now. Since every DFPT U hook (cal_docc, build_dv_u,
    // dftu_onsite, born/docc contractions) is a no-op reservation, running
    // anyway would converge cleanly while silently dropping the whole
    // first-order U response; init() must therefore reject the run
    // explicitly (same fail-loud pattern as the metallic-sampling guard).
    // The accessor semantics (with_u true, u_active following the
    // occupation-matrix state) are pinned separately in
    // DFPT_PW_DataTest.DftuReservationProviderUsability.
    Plus_U_Base dftu;
    dfpt.set_qmesh(1, 1, 1);
    psi::Psi<std::complex<double>> psi;
    // death tests match the child's stderr, while WARNING_QUIT writes the
    // NOTICE block to std::cout; bridge the two inside the statement
    EXPECT_EXIT({
        std::cout.rdbuf(std::cerr.rdbuf());
        dfpt.init(ucell, psi, nullptr, nullptr, nullptr, std::vector<double>(),
                  ModuleBase::matrix(), ModuleBase::matrix(), nullptr, 1.0, 15.0, &dftu);
    }, ::testing::ExitedWithCode(1), "DFT\\+U with DFPT is not supported");
}
