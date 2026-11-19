#include "source_cell/dom_decomp.h"
#include "source_base/parallel_cell.h"
#include "gtest/gtest.h"
#include "source_md/unittests/setcell.h"
#include "lj_test_solver.h"
#include "source_io/module_parameter/parameter.h"
#include "source_md/md_func.h"
#include <fstream>
#include <cstdint>
#include <string>
constexpr double doublethreshold = 1e-12;

// LJ force-field parameters and numerical results, including the MDCell path.

class LJ_pot_test : public testing::Test
{
  protected:
    ModuleBase::matrix stress;
    double potential;
    int natom;
    UnitCell ucell;
    Input_para inp;
    std::ofstream running_log;
    std::string output_dir;

    void SetUp()
    {
        Setcell::setupcell(ucell);

        natom = ucell.nat;
        stress.create(3, 3);

        Setcell::parameters(inp);
        output_dir = testing::TempDir();
    }

    void TearDown()
    {
    }
};

TEST_F(LJ_pot_test, potential)
{
    ModuleESolver::ESolver* p_esolver = make_lj_test_solver();
    MDCell mdcell;
    mdcell = Setcell::setup_mdcell(ucell);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), mdcell.latvec(), mdcell.lat0(), 0.0, 0.0);
    p_esolver->before_all_runners(mdcell, inp);
    const bool positions_updated = false;
    MD_func::force_virial(p_esolver, 0, positions_updated, mdcell, decomp, potential, true, stress);
    EXPECT_NEAR(potential, -0.011957818623534381, doublethreshold);
}

TEST_F(LJ_pot_test, unitcell_compatibility)
{
    ModuleESolver::ESolver* p_esolver = make_lj_test_solver();
    p_esolver->before_all_runners(ucell, inp);
    p_esolver->runner(ucell, 0);
    p_esolver->cal_force(ucell, stress);

    EXPECT_NEAR(p_esolver->cal_energy(), -0.023915637247068761, doublethreshold);
    EXPECT_NEAR(stress(0, 0), 0.00099635466178755387, doublethreshold);

    delete p_esolver;
}

TEST_F(LJ_pot_test, force)
{
    ModuleESolver::ESolver* p_esolver = make_lj_test_solver();
    MDCell mdcell;
    mdcell = Setcell::setup_mdcell(ucell);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), mdcell.latvec(), mdcell.lat0(), 0.0, 0.0);
    p_esolver->before_all_runners(mdcell, inp);
    const bool positions_updated = false;
    MD_func::force_virial(p_esolver, 0, positions_updated, mdcell, decomp, potential, true, stress);
    const std::vector<LocalAtom>& atoms = mdcell.owned_atoms();
    EXPECT_NEAR(atoms[0].force.x, 0.00049817733089377704, doublethreshold);
    EXPECT_NEAR(atoms[0].force.y, 0.00082237246837022328, doublethreshold);
    EXPECT_NEAR(atoms[0].force.z, 0.0, doublethreshold);
    EXPECT_NEAR(atoms[1].force.x, -0.00064758615201580339, doublethreshold);
    EXPECT_NEAR(atoms[1].force.y, -0.00066924999462089304, doublethreshold);
    EXPECT_NEAR(atoms[1].force.z, 0.0, doublethreshold);
    EXPECT_NEAR(atoms[2].force.x, -0.00035411224839165616, doublethreshold);
    EXPECT_NEAR(atoms[2].force.y, 0.0008091080910885112, doublethreshold);
    EXPECT_NEAR(atoms[2].force.z, 0.0, doublethreshold);
    EXPECT_NEAR(atoms[3].force.x, 0.00050352106951368229, doublethreshold);
    EXPECT_NEAR(atoms[3].force.y, -0.00096223056483784122, doublethreshold);
    EXPECT_NEAR(atoms[3].force.z, 0.0, doublethreshold);
}

TEST_F(LJ_pot_test, mdcell_cal_force)
{
    std::unique_ptr<ModuleESolver::ESolver> p_esolver(make_lj_test_solver());
    MDCell mdcell;
    mdcell = Setcell::setup_mdcell(ucell);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), mdcell.latvec(), mdcell.lat0(), 0.0, 0.0);
    p_esolver->before_all_runners(mdcell, inp);
    decomp.prepare_neighbors(mdcell);
    p_esolver->runner(mdcell, 0);

    ModuleBase::matrix force;
    p_esolver->cal_force(mdcell, force);
    for (int iat = 0; iat < mdcell.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(iat)];
        EXPECT_DOUBLE_EQ(force(iat, 0), atom.force.x);
        EXPECT_DOUBLE_EQ(force(iat, 1), atom.force.y);
        EXPECT_DOUBLE_EQ(force(iat, 2), atom.force.z);
    }
}

TEST_F(LJ_pot_test, stress)
{
    ModuleESolver::ESolver* p_esolver = make_lj_test_solver();
    MDCell mdcell;
    mdcell = Setcell::setup_mdcell(ucell);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), mdcell.latvec(), mdcell.lat0(), 0.0, 0.0);
    p_esolver->before_all_runners(mdcell, inp);
    const bool positions_updated = false;
    MD_func::force_virial(p_esolver, 0, positions_updated, mdcell, decomp, potential, true, stress);
    EXPECT_NEAR(stress(0, 0), 8.0360222227631859e-07, doublethreshold);
    EXPECT_NEAR(stress(0, 1), 1.7207745586539077e-07, doublethreshold);
    EXPECT_NEAR(stress(0, 2), 0.0, doublethreshold);
    EXPECT_NEAR(stress(1, 0), 1.7207745586539077e-07, doublethreshold);
    EXPECT_NEAR(stress(1, 1), 1.0630708613186662e-06, doublethreshold);
    EXPECT_NEAR(stress(1, 2), 0.0, doublethreshold);
    EXPECT_NEAR(stress(2, 0), 0.0, doublethreshold);
    EXPECT_NEAR(stress(2, 1), 0.0, doublethreshold);
    EXPECT_NEAR(stress(2, 2), 6.4275429572682057e-07, doublethreshold);
}

TEST_F(LJ_pot_test, RcutSearchRadius)
{
    ModuleMD::LJ* p_esolver = new ModuleMD::LJ(output_dir, running_log);
    ucell.ntype = 2;
    std::vector<double> rcut = {3.0};
    p_esolver->rcut_search_radius_for_testing(ucell.ntype, rcut);

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j < ucell.ntype; j++)
        {
            EXPECT_NEAR(p_esolver->get_lj_rcut()(i, j), 3.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
        }
    }
    EXPECT_NEAR(p_esolver->get_search_radius(), 3.0 * ModuleBase::ANGSTROM_AU + 0.01, doublethreshold);

    rcut = {3.0, 4.0, 5.0};
    p_esolver->rcut_search_radius_for_testing(ucell.ntype, rcut);
    EXPECT_NEAR(p_esolver->get_lj_rcut()(0, 0), 3.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->get_lj_rcut()(0, 1), 4.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->get_lj_rcut()(1, 0), 4.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->get_lj_rcut()(1, 1), 5.0 * ModuleBase::ANGSTROM_AU, doublethreshold);
    EXPECT_NEAR(p_esolver->get_search_radius(), 5.0 * ModuleBase::ANGSTROM_AU + 0.01, doublethreshold);
}

TEST_F(LJ_pot_test, SetC6C12)
{
    ModuleMD::LJ* p_esolver = new ModuleMD::LJ(output_dir, running_log);
    ucell.ntype = 2;

    // no rule
    int rule = 1;
    std::vector<double> lj_epsilon = {0.1, 0.2, 0.3};
    std::vector<double> lj_sigma = {0.2, 0.4, 0.6};

    p_esolver->set_c6_c12_for_testing(ucell.ntype, rule, lj_epsilon, lj_sigma);

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            int k = i * (i + 1) / 2 + j;
            double temp = pow(lj_sigma[k] * ModuleBase::ANGSTROM_AU, 6);
            EXPECT_NEAR(p_esolver->get_lj_c6()(i, j), 4.0 * lj_epsilon[k] * temp / ModuleBase::Ry_to_eV, doublethreshold);
            EXPECT_NEAR(p_esolver->get_lj_c12()(i, j), p_esolver->get_lj_c6()(i, j) * temp, doublethreshold);
            EXPECT_DOUBLE_EQ(p_esolver->get_lj_c6()(i, j), p_esolver->get_lj_c6()(j, i));
            EXPECT_DOUBLE_EQ(p_esolver->get_lj_c12()(i, j), p_esolver->get_lj_c12()(j, i));
        }
    }

    // rule 1
    rule = 1;
    lj_epsilon = {0.1, 0.2};
    lj_sigma = {0.2, 0.4};

    p_esolver->set_c6_c12_for_testing(ucell.ntype, rule, lj_epsilon, lj_sigma);

    for (int i = 0; i < ucell.ntype; i++)
    {
        double temp = pow(lj_sigma[i] * ModuleBase::ANGSTROM_AU, 6);
        EXPECT_NEAR(p_esolver->get_lj_c6()(i, i), 4.0 * lj_epsilon[i] * temp / ModuleBase::Ry_to_eV, doublethreshold);
        EXPECT_NEAR(p_esolver->get_lj_c12()(i, i), p_esolver->get_lj_c6()(i, i) * temp, doublethreshold);

        for (int j = 0; j < i; j++)
        {
            EXPECT_NEAR(p_esolver->get_lj_c6()(i, j),
                        std::sqrt(p_esolver->get_lj_c6()(i, i) * p_esolver->get_lj_c6()(j, j)),
                        doublethreshold);
            EXPECT_NEAR(p_esolver->get_lj_c12()(i, j),
                        std::sqrt(p_esolver->get_lj_c12()(i, i) * p_esolver->get_lj_c12()(j, j)),
                        doublethreshold);
            EXPECT_DOUBLE_EQ(p_esolver->get_lj_c6()(i, j), p_esolver->get_lj_c6()(j, i));
            EXPECT_DOUBLE_EQ(p_esolver->get_lj_c12()(i, j), p_esolver->get_lj_c12()(j, i));
        }
    }

    // rule 2
    rule = 2;
    lj_epsilon = {0.1, 0.2};
    lj_sigma = {0.2, 0.4};

    p_esolver->set_c6_c12_for_testing(ucell.ntype, rule, lj_epsilon, lj_sigma);

    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double temp = pow((lj_sigma[i] + lj_sigma[j]) / 2 * ModuleBase::ANGSTROM_AU, 6);
            EXPECT_NEAR(p_esolver->get_lj_c6()(i, j),
                        4.0 * std::sqrt(lj_epsilon[i] * lj_epsilon[j]) * temp / ModuleBase::Ry_to_eV,
                        doublethreshold);
            EXPECT_NEAR(p_esolver->get_lj_c12()(i, j), p_esolver->get_lj_c6()(i, j) * temp, doublethreshold);
            EXPECT_DOUBLE_EQ(p_esolver->get_lj_c6()(i, j), p_esolver->get_lj_c6()(j, i));
            EXPECT_DOUBLE_EQ(p_esolver->get_lj_c12()(i, j), p_esolver->get_lj_c12()(j, i));
        }
    }
}

TEST_F(LJ_pot_test, CalEnShift)
{
    ModuleMD::LJ* p_esolver = new ModuleMD::LJ(output_dir, running_log);
    ucell.ntype = 2;

    std::vector<double> rcut = {3.0};
    p_esolver->rcut_search_radius_for_testing(ucell.ntype, rcut);

    int rule = 1;
    std::vector<double> lj_epsilon = {0.1, 0.2, 0.3};
    std::vector<double> lj_sigma = {0.2, 0.4, 0.6};
    p_esolver->set_c6_c12_for_testing(ucell.ntype, rule, lj_epsilon, lj_sigma);

    // false
    p_esolver->cal_en_shift_for_testing(ucell.ntype, false);
    for (int i = 0; i < ucell.ntype; i++)
    {
        for (int j = 0; j < ucell.ntype; j++)
        {
            EXPECT_DOUBLE_EQ(p_esolver->get_en_shift()(i, j), 0.0);
        }
    }

    // true
    p_esolver->cal_en_shift_for_testing(ucell.ntype, true);
    EXPECT_NEAR(p_esolver->get_en_shift()(0, 0), -2.5810212013100967e-09, doublethreshold);
    EXPECT_NEAR(p_esolver->get_en_shift()(0, 1), -3.303688865319793e-07, doublethreshold);
    EXPECT_NEAR(p_esolver->get_en_shift()(1, 0), -3.303688865319793e-07, doublethreshold);
    EXPECT_NEAR(p_esolver->get_en_shift()(1, 1), -5.6443326024140752e-06, doublethreshold);
}

#include "source_md/module_lj/lj_forcefield.h"
#include "source_esolver/esolver_forcefields.h"
#include "source_cell/dom_decomp.h"
#include "source_md/unittests/setcell.h"

#include "gtest/gtest.h"

#include <memory>
#include <fstream>
#include <utility>

TEST(LJPotential, UnitCellEnergyForceAndStress)
{
    UnitCell cell;
    Setcell::setupcell(cell);
    Input_para input;
    Setcell::parameters(input);
    const std::string output_dir = testing::TempDir();
    std::ofstream running_log;
    std::unique_ptr<ModuleMD::ForceField> potential(new ModuleMD::LJ(output_dir, running_log));
    ModuleESolver::esolver_forcefields solver(std::move(potential), running_log);
    solver.before_all_runners(cell, input);
    solver.runner(cell, 0);

    EXPECT_NEAR(solver.cal_energy(), -0.023915637247068761, 1e-12);
    ModuleBase::matrix force;
    solver.cal_force(cell, force);
    EXPECT_NEAR(force(0, 0), 0.00099635466178755387, 1e-12);
    ModuleBase::matrix stress;
    solver.cal_stress(cell, stress);
    EXPECT_NEAR(stress(0, 0), 1.6072044445526372e-6, 1e-10);
}

TEST(LJPotential, MDCellWritesOwnedForces)
{
    UnitCell unitcell;
    Setcell::setupcell(unitcell);
    MDCell cell = Setcell::setup_mdcell(unitcell);
    Input_para input;
    Setcell::parameters(input);
    const std::string output_dir = testing::TempDir();
    std::ofstream running_log;
    std::unique_ptr<ModuleMD::ForceField> potential(new ModuleMD::LJ(output_dir, running_log));
    ModuleESolver::esolver_forcefields solver(std::move(potential), running_log);
    solver.before_all_runners(cell, input);
    DomainDecomposition decomp;
    decomp.init(ModuleBase::world_comm_domain(), cell.latvec(), cell.lat0(), 0.0, 0.0);
    decomp.prepare_neighbors(cell);
    solver.runner(cell, 0);

    EXPECT_NEAR(solver.cal_energy(), -0.023915637247068761, 1e-12);
    ModuleBase::matrix force;
    solver.cal_force(cell, force);
    for (std::size_t i = 0; i < cell.owned_atoms().size(); ++i)
    {
        EXPECT_DOUBLE_EQ(force(i, 0), cell.owned_atoms()[i].force.x);
        EXPECT_DOUBLE_EQ(force(i, 1), cell.owned_atoms()[i].force.y);
        EXPECT_DOUBLE_EQ(force(i, 2), cell.owned_atoms()[i].force.z);
    }
}

#ifdef __MPI
TEST(LJPotential, SplitCommunicatorReducesWithinCell)
{
    UnitCell unitcell;
    Setcell::setupcell(unitcell);
    Input_para input;
    Setcell::parameters(input);

    ModuleBase::CommunicationDomain domain;
    domain.initialize(MPI_COMM_SELF);
    std::vector<LocalAtom> atoms;
    for (int index = 0; index < unitcell.atoms[0].na; ++index)
    {
        const ModuleBase::Vector3<double> zero_force(0.0, 0.0, 0.0);
        const double mass_au = unitcell.atoms[0].mass / ModuleBase::AU_to_MASS;
        const LocalAtom atom(unitcell.atoms[0].tau[index],
                             unitcell.atoms[0].taud[index],
                             unitcell.atoms[0].vel[index],
                             zero_force,
                             unitcell.atoms[0].mbl[index],
                             mass_au,
                             0,
                             index,
                             0);
        atoms.push_back(atom);
    }
    const std::vector<std::string> labels(1, unitcell.atoms[0].label);
    const std::vector<double> masses(1, unitcell.atoms[0].mass);
    const std::vector<std::int64_t> counts(1, unitcell.atoms[0].na);
    MDCell cell;
    cell.initialize_from_owned_atoms(unitcell.latvec,
                                     unitcell.GT,
                                     unitcell.lat0,
                                     unitcell.omega,
                                     unitcell.nat,
                                     atoms,
                                     labels,
                                     masses,
                                     counts,
                                     0.0,
                                     domain);

    std::ofstream running_log;
    const std::string output_dir = testing::TempDir();
    ModuleMD::LJ potential(output_dir, running_log);
    potential.initialize(cell, input);
    DomainDecomposition decomposition;
    decomposition.init(domain, cell.latvec(), cell.lat0(), 0.0, 0.0);
    decomposition.prepare_neighbors(cell);
    double energy = 0.0;
    ModuleBase::matrix force;
    ModuleBase::matrix stress(3, 3);
    potential.compute(cell, energy, force, stress);

    EXPECT_NEAR(energy, -0.023915637247068761, 1e-12);
    EXPECT_NEAR(stress(0, 0), 1.6072044445526372e-6, 1e-10);
}
#endif
