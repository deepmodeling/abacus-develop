#include <fstream>
#include <iterator>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_cell/unitcell.h"
#include "source_relax/ions_move_basic.h"
#include "source_relax/ions_move_lbfgs.h"
#include "source_relax/relax_criteria.h"

/************************************************
 *  unit tests of class Ions_Move_LBFGS
 ***********************************************/

/**
 * - Tested Functions:
 *   - Ions_Move_LBFGS::relax_step() converged / not-converged branches
 *   - Ions_Move_LBFGS::relax_step() fills Ions_Move_Basic::largest_grad
 *     (the value reported by IonCellOptimizer::relax_step in the running log)
 */

namespace unitcell
{
// Stub: the real update_pos_tau performs lattice-vector algebra that needs a
// fully built UnitCell; these tests only care about the convergence decision.
void update_pos_tau(const Lattice&, const double*, const int, const int, Atom*)
{
}
} // namespace unitcell

UnitCell::UnitCell()
{
    Coordinate = "Direct";
    latName = "none";
    lat0 = 1.0;
    latvec.Identity();

    ntype = 1;
    nat = 1;
    iat2it = nullptr;
    iat2ia = nullptr;
    iwt2iat = nullptr;
    iwt2iw = nullptr;
    itia2iat.create(1, 1);

    atoms = new Atom[ntype];
    set_atom_flag = true;
    atoms[0].label = "Si";
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
Atom::Atom()
{
    na = 1;
    tau.resize(na);
    dis.resize(na);
    mag.resize(na);
    mbl.resize(na);
    vel.resize(na);
    taud.resize(na);
    mbl[0] = {1, 1, 1};
}
Atom::~Atom()
{
}
#ifdef __MPI
void Atom::bcast_atom()
{
}
void Atom::bcast_atom2()
{
}
#endif
Atom_pseudo::Atom_pseudo()
{
}
Atom_pseudo::~Atom_pseudo()
{
}
pseudo::pseudo()
{
}
pseudo::~pseudo()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}
int ModuleSymmetry::Symmetry::symm_flag = 0;
void ModuleSymmetry::Symmetry::symmetrize_mat3(ModuleBase::matrix& sigma, const Lattice& lat) const
{
}
void ModuleSymmetry::Symmetry::symmetrize_vec3_nat(double* v) const
{
}

// Friend of Parameter; the only sanctioned way to write PARAM.input in tests.
// Not needed anymore: allocate() takes relax_bfgs_rmax and out_level explicitly.
class IonsMoveLBFGSTest : public ::testing::Test
{
  protected:
    Ions_Move_LBFGS lbfgs;
    const int natom = 1;

    void SetUp() override
    {
        lbfgs.allocate(natom, 0.2, "m"); // keep stdout quiet
    }
};

// A force below force_thr_ev must be reported as converged, with
// Ions_Move_Basic::largest_grad filled for the running-log report.
TEST_F(IonsMoveLBFGSTest, RelaxStepConverged)
{
    UnitCell ucell;
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 1.0e-4; // Ry/Bohr, ~0.005 eV/Angstrom
    const double etot = 0.0;
    std::ofstream ofs("lbfgs_converged.log");
    Relax_Criteria criteria;
    criteria.force_thr_ev = 1.0; // eV/Angstrom

    const bool converged = lbfgs.relax_step(force, ucell, etot, ofs, criteria);
    ofs.close();
    std::remove("lbfgs_converged.log");

    EXPECT_TRUE(converged);
    EXPECT_NEAR(Ions_Move_Basic::largest_grad, 1.0e-4, 1e-12); // Ry/Bohr, lat0 = 1
    EXPECT_LT(Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A, criteria.force_thr_ev);
}

// A force above force_thr_ev must be reported as not converged.
TEST_F(IonsMoveLBFGSTest, RelaxStepNotConverged)
{
    UnitCell ucell;
    ModuleBase::matrix force(natom, 3);
    force(0, 0) = 0.1; // Ry/Bohr, ~51.4 eV/Angstrom
    const double etot = 0.0;
    std::ofstream ofs("lbfgs_not_converged.log");
    Relax_Criteria criteria;
    criteria.force_thr_ev = 1.0; // eV/Angstrom

    const bool converged = lbfgs.relax_step(force, ucell, etot, ofs, criteria);
    ofs.close();
    std::remove("lbfgs_not_converged.log");

    EXPECT_FALSE(converged);
    EXPECT_NEAR(Ions_Move_Basic::largest_grad, 0.1, 1e-12);
    EXPECT_GT(Ions_Move_Basic::largest_grad * ModuleBase::Ry_to_eV / ModuleBase::BOHR_TO_A, criteria.force_thr_ev);
}
