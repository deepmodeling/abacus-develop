#include "source_lcao/module_ri/sternheimer_abacus_fd_adapter.h"

#include "source_base/matrix3.h"

#include <complex>
#include <gtest/gtest.h>
#include <memory>
#include <stdexcept>
#include <vector>

TEST(SternheimerABACUSFDAdapter, BuildsOrthogonalGridFromABACUSLattice)
{
    const ModuleBase::Matrix3 latvec(3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0);
    const auto grid_data = ModuleRI::make_sternheimer_fd_grid_from_lattice(6, 8, 10, 6 * 8 * 10, 2.0, latvec);

    EXPECT_EQ(grid_data.grid.nx, 6);
    EXPECT_EQ(grid_data.grid.ny, 8);
    EXPECT_EQ(grid_data.grid.nz, 10);
    EXPECT_TRUE(grid_data.grid.periodic);
    EXPECT_DOUBLE_EQ(grid_data.grid.hx, 1.0);
    EXPECT_DOUBLE_EQ(grid_data.grid.hy, 1.0);
    EXPECT_DOUBLE_EQ(grid_data.grid.hz, 1.0);
    EXPECT_DOUBLE_EQ(grid_data.volume_element, 1.0);
}

TEST(SternheimerABACUSFDAdapter, RejectsDistributedRealSpaceGridForDensePrototype)
{
    const ModuleBase::Matrix3 latvec(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    EXPECT_THROW(ModuleRI::make_sternheimer_fd_grid_from_lattice(4, 4, 4, 32, 8.0, latvec), std::invalid_argument);
}

TEST(SternheimerABACUSFDAdapter, RejectsNonOrthogonalLatticeForCurrentStencil)
{
    const ModuleBase::Matrix3 latvec(1.0, 0.2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    EXPECT_THROW(ModuleRI::make_sternheimer_fd_grid_from_lattice(4, 4, 4, 64, 8.0, latvec), std::invalid_argument);
}

TEST(SternheimerABACUSFDAdapter, BuildsHamiltonianWithNonlocalProjector)
{
    ModuleRI::SternheimerABACUSFDGridData grid_data;
    grid_data.grid.nx = 2;
    grid_data.grid.ny = 1;
    grid_data.grid.nz = 1;
    grid_data.grid.hx = 1.0;
    grid_data.grid.hy = 1.0;
    grid_data.grid.hz = 1.0;
    grid_data.grid.periodic = true;
    grid_data.volume_element = 1.0;

    using Projector = ModuleRI::SternheimerFDNonlocalProjector;
    using Complex = Projector::Complex;
    Projector::ProjectorBlock block;
    block.projectors = {{Complex(1.0, 0.0), Complex(0.0, 0.0)}};
    block.d_matrix = {{Complex(2.0, 0.0)}};
    const auto nonlocal_projector = std::make_shared<Projector>(2, grid_data.volume_element, std::vector<Projector::ProjectorBlock>{block});

    const auto hamiltonian = ModuleRI::make_sternheimer_fd_hamiltonian_from_local_potential(
        grid_data,
        std::vector<double>{0.0, 0.0},
        0.0,
        nonlocal_projector);

    ASSERT_NE(hamiltonian.nonlocal_projector(), nullptr);

    const ModuleRI::SternheimerFDHamiltonian::Vector psi = {Complex(3.0, 0.0), Complex(4.0, 0.0)};
    ModuleRI::SternheimerFDHamiltonian::Vector hpsi;
    hamiltonian.apply(psi, hpsi);

    ASSERT_EQ(hpsi.size(), psi.size());
    EXPECT_NEAR(hpsi[0].real(), 6.0, 1.0e-12);
    EXPECT_NEAR(hpsi[0].imag(), 0.0, 1.0e-12);
    EXPECT_NEAR(hpsi[1].real(), 0.0, 1.0e-12);
    EXPECT_NEAR(hpsi[1].imag(), 0.0, 1.0e-12);
}
