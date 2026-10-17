#include <gtest/gtest.h>
#include "source_lcao/module_lr/state_track.h"
#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include "source_psi/psi.h"
#include <memory>
#include <sstream>

TEST(StateTrack, ClosedShellTripletIgnoresDownSpinVirtualRotation)
{
    // Production follow_cross_root with an identity AO overlap from the existing
    // two-center mock: well-separated one-orbital atoms have no interatom pairs.
    UnitCell cell;
    std::unique_ptr<Atom> atoms(new Atom);
    cell.atoms = atoms.get();
    cell.ntype = 1;
    cell.nat = 3;
    cell.lat0 = 1.0;
    cell.latvec = ModuleBase::Matrix3(10, 0, 0, 0, 10, 0, 0, 0, 10);
    cell.iat2it = {0, 0, 0};
    cell.iat2ia = {0, 1, 2};
    atoms->nw = 1;
    atoms->na = 3;
    atoms->iw2l = {0};
    atoms->iw2n = {0};
    atoms->iw2m = {0};
    atoms->tau = {{0, 0, 0}, {2, 0, 0}, {4, 0, 0}};
    const std::vector<int> offsets{0, 1, 2};
    cell.set_iat2iwt_for_test(offsets, 1);
    Parallel_2D pc;
#ifdef __MPI
    pc.init(3, 3, 1, MPI_COMM_WORLD);
#else
    pc.set_serial(3, 3);
#endif
    Parallel_Orbitals pmat;
    std::vector<Parallel_2D> px(2);
#ifdef __MPI
    pmat.set(3, 3, 1, pc.blacs_ctxt);
    for (auto& block: px)
    {
        block.set(2, 1, 1, pc.blacs_ctxt);
    }
#else
    pmat.set_serial(3, 3);
    for (auto& block: px)
    {
        block.set_serial(2, 1);
    }
#endif
    pmat.set_atomic_trace(cell.get_iat2iwt(), cell.nat, 3);
    const std::vector<int> ngk{3, 3};
    psi::Psi<double> orbitals(2, 3, 3, ngk, true);
    orbitals.zero_out();
    for (int spin = 0; spin < 2; ++spin)
        for (int band = 0; band < 3; ++band)
        {
            orbitals(spin, band, band) = 1.0;
        }
    TwoCenterIntegrator integrator;
    const std::vector<double> cutoff{0.5};
    const std::vector<int> nocc{1, 1};
    const std::vector<int> nvirt{2, 2};
    const LR::RootInputs<double> inputs{cell, cutoff, integrator, orbitals, pc, pmat, px,
        nocc, nvirt, false, 1, false};
    const std::vector<double> amplitudes{1, 0, 0, 1};
    const std::vector<double> energies{1, 2};
    LR::RootBasis<double> basis;
    std::vector<double> previous;
    std::ostringstream log;
    int target = 0;
    LR::follow_cross_root(inputs, amplitudes.data(), 2, 2, 0, energies, 0.0, target, previous, basis, log);
    // Rotate only spin-down orbitals; closed-shell triplet X and spin-up C stay fixed.
    orbitals(1, 1, 1) = 0.0;
    orbitals(1, 2, 2) = 0.0;
    orbitals(1, 1, 2) = 1.0;
    orbitals(1, 2, 1) = 1.0;
    LR::follow_cross_root(inputs, amplitudes.data(), 2, 2, 0, energies, 0.0, target, previous, basis, log);
    EXPECT_EQ(target, 0);
    ASSERT_EQ(previous.size(), 2u);
    EXPECT_NEAR(previous[0], 1.0, 1e-13);
    EXPECT_NEAR(previous[1], 0.0, 1e-13);
}
