#ifndef DFPT_SERIAL_FIXTURE_H
#define DFPT_SERIAL_FIXTURE_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/qlist.h"
#include "source_cell/unitcell.h"
#include "source_pw/module_dfpt/dfpt_pw_data.h"

#include "gtest/gtest.h"
#include <complex>

// Shared serial-side gtest fixture for the DFPT unit tests
// (dfpt_pert/rho/phon/q0_serial_test.cpp). Everything runs without
// __MPI: the plane-wave bases are built through the real serial
// initgrids/initparameters/setuptransform path on a shared FFT grid,
// exactly like the production setup_pwrho/setup_pwwfc sequence.
//
// All members the fixture touches (UnitCell geometry fields, the
// Atom/pseudo public data, QList::nkstot / kvec_d and the
// DFPT_PW_Data::init entry) are public, so the tests include the
// cell/qlist/dfpt headers normally.

class DFPTSerialBase : public testing::Test
{
  protected:
    const double lat0_ = 1.8897261254578281;
    const double ecutwfc_ = 2.5; // Ry
    // rho cutoff inflated to 9x ecutwfc so every Delta = G''-G' of the
    // convolution lies inside the rho ball and nothing aliases
    const double rho_mult_ = 9.0;
    const double a_ = 10.0; // cubic edge in lat0 units

    ModuleBase::Matrix3 latvec_;
    UnitCell ucell_;
    ModulePW::PW_Basis pw_rho_;
    ModulePW::PW_Basis_K pw_wfc_;
    ModuleCell::QList qlist_;
    ModuleDFPT::DFPT_PW_Data data_;
    ModuleBase::Matrix3 G_; // = latvec^-T, the reciprocal builder

    // default (k, q) of the pert/phon/rho fixtures: q is generic and
    // k = -q so k+q = 0: the k+q ball then stays inside the ground-state
    // G list (single-k limitation documented in DFPT_KQ_Basis)
    const ModuleBase::Vector3<double> q_d_{0.13, 0.0, 0.07};
    const ModuleBase::Vector3<double> k_d_{-0.13, 0.0, -0.07};
    ModuleBase::Vector3<double> q_cart_;
    const ModuleBase::Vector3<double> tau_{1.1, 2.3, 0.7}; // lat0 units

    // single-atom Coulomb cell + shared-grid bases at the default (k_d_, q_d_)
    void SetUp() override;
    void TearDown() override;

    // cubic single-atom Coulomb cell (iat2it/iat2ia allocated)
    void SetUpCell();

    // (re)initialize the bases and the shared data wiring for a given
    // (k, q) pair and band count; SetUp uses the default fixture values
    void SetupBases(const ModuleBase::Vector3<double>& k_d, const ModuleBase::Vector3<double>& q_d, int nbands);

    void MakeCoulombAtom();
    void MakeNCAtom();

    /// extend the two-atom cell of MakeTwoAtomCell with a nonlinear core
    /// correction on the second species: a Gaussian rho_atc of total charge
    /// nlcc_core_charge_ on the existing radial grid (msh points). The
    /// atom sits at the non-origin tau_nlcc_ = (0.25, 0.31, 0.17) -- the
    /// anti-blind-spot position: every structure-factor phase is active,
    /// so a mirrored phase sign (exp(+iG.tau) instead of the QE addcore
    /// exp(-iG.tau)) places the response at the wrong site and is caught
    /// directly by the drhoc-localization test.
    void MakeNlccTwoAtomCell();

    /// total integrated core charge (in e) of the nlcc Gaussian above
    static constexpr double nlcc_core_charge_ = 0.24;
    const ModuleBase::Vector3<double> tau_nlcc_{0.25, 0.31, 0.17}; // lat0 units
    double nlcc_sigma_ = 1.5; // bohr, grid-resolvable (spacing ~0.63 bohr)

    // reconfigure the cell as a two-atom Z=4/Z=2 crystal breaking all symmetry
    void MakeTwoAtomCell();

    // key of an integer FFT triple (gcar * a is integral on the cubic cell)
    long long FKey(int ix, int iy, int iz) const;
    long long GKey(const ModuleBase::Vector3<double>& g) const;

    // analytic Coulomb local potential (Ry) at |g|^2 in bohr^-2, mirroring
    // vl_pw.cpp::vloc_coulomb independently of DFPT_Pert::vloc_at_g
    double VlocCoulomb(double g2_bohr) const;

    // analytic dVloc/dtau_alpha coefficient at displacement vector w (1/lat0);
    // GS structure-factor convention (stru_fac.cpp): exp(-i 2pi (g.tau)) and
    // dV/dtau = -i (Delta+q)_alpha tpiba Vloc exp(-i 2pi (Delta+q).tau)
    std::complex<double> AnalyticDVloc(int dir, const ModuleBase::Vector3<double>& w) const;
};

#endif // DFPT_SERIAL_FIXTURE_H
