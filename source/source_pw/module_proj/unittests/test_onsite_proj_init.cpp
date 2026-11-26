#include "source_pw/module_proj/onsite_proj.h"
#include "source_pw/module_proj/radial_proj.h"

#include <gtest/gtest.h>

#include <cmath>

namespace
{
int tools_created = 0;
int tools_destroyed = 0;
int tools_npw = 0;
const void* tools_psi = nullptr;
const void* tools_buffer = nullptr;
}

// Observe the force/stress backend lifetime and borrowed inputs without running
// FFTs or force kernels. The projector's init, SBT and buffer allocation are real.
namespace hamilt
{
template <typename T, typename Device>
Onsite_Proj_tools<T, Device>::Onsite_Proj_tools(
    const std::vector<int>&,
    const std::vector<int>&,
    const ModuleBase::realArray&,
    const ModuleBase::matrix&,
    std::complex<T>* buffer,
    const UnitCell*,
    const psi::Psi<std::complex<T>, Device>* wavefunction,
    const K_Vectors*,
    const ModulePW::PW_Basis_K* basis,
    const Structure_Factor*,
    const ModuleBase::matrix&,
    const ModuleBase::matrix&)
{
    ++tools_created;
    tools_npw = basis->npwk_max;
    tools_psi = wavefunction;
    tools_buffer = buffer;
}

template <typename T, typename Device>
Onsite_Proj_tools<T, Device>::~Onsite_Proj_tools()
{
    ++tools_destroyed;
}

template <typename T, typename Device>
void Onsite_Proj_tools<T, Device>::refresh_table()
{
}

template class Onsite_Proj_tools<double, base_device::DEVICE_CPU>;
#if defined(__CUDA) || defined(__ROCM)
template class Onsite_Proj_tools<double, base_device::DEVICE_GPU>;
#endif
}

// Structure factors are only borrowed during init; no geometry evaluation is needed.
Structure_Factor::Structure_Factor() = default;
Structure_Factor::~Structure_Factor() = default;

namespace projectors
{
class OnsiteProjectorInitTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        tools_created = 0;
        tools_destroyed = 0;
        cell.omega = 8.0;
        basis.npwk_max = 4;
        kv.isk = {0};
        // Seed the state after orbital loading so no pseudopotential/orbital
        // files or SCF calculation are needed to exercise reinitialization.
        projector.initialed = true;
        projector.ucell = &cell;
        projector.ntype = 1;
        projector.tot_nproj = 1;
        projector.it2iproj = {{0}};
        projector.lproj = {0};
        projector.rgrid = {0.0, 0.25, 0.5, 0.75, 1.0};
        projector.projs = {{1.0, 0.8, 0.5, 0.2, 0.0}};
        RadialProjection::build_sbt_tab({1}, projector.rgrid, projector.projs,
                                      projector.lproj, 12, 0.1, cell.omega, 1,
                                      projector.tab, projector.nhtol);
        projector.npwx_ = basis.npwk_max;
        projector.size_vproj = basis.npwk_max;
        projector.isk_ = kv.isk.data();
        base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_CPU>()(
            projector.tab_atomic_, projector.size_vproj, "OnsiteP::tab_atomic_");
        projector.fs_tools = new hamilt::Onsite_Proj_tools<double, base_device::DEVICE_CPU>(
            {1}, projector.lproj, projector.tab, projector.nhtol, projector.tab_atomic_,
            &cell, &wavefunction, &kv, &basis, &sf, wg, ekb);
        cell.cell_parameter_updated = false;
    }

    void initialize(const psi::Psi<std::complex<double>>& psi)
    {
        projector.init("", &cell, psi, kv, basis, sf, 0.5, 8, 0.1, wg, ekb);
    }

    double radial_zero() const
    {
        return projector.tab(0, 0, 0);
    }

    void check_state(int capacity, int plane_waves)
    {
        EXPECT_EQ(projector.tab.getBound3(), capacity);
        EXPECT_EQ(projector.get_npwx(), plane_waves);
        EXPECT_EQ(projector.size_vproj, plane_waves);
        ASSERT_NE(projector.get_tab_atomic(), nullptr);
        EXPECT_EQ(tools_npw, plane_waves);
        EXPECT_EQ(tools_buffer, projector.get_tab_atomic());
        EXPECT_EQ(projector.isk_, kv.isk.data());
    }

    UnitCell cell;
    ModulePW::PW_Basis_K basis;
    K_Vectors kv;
    Structure_Factor sf;
    psi::Psi<std::complex<double>> wavefunction;
    ModuleBase::matrix wg;
    ModuleBase::matrix ekb;
    OnsiteProjector<double, base_device::DEVICE_CPU> projector;
};

TEST_F(OnsiteProjectorInitTest, ChangedVolumeRebuildsTableAndRebindsTools)
{
    ASSERT_EQ(tools_created, 1);
    const double initial = radial_zero();
    ASSERT_GT(std::abs(initial), 0.0);
    psi::Psi<std::complex<double>> next_wavefunction;
    cell.omega = 32.0;
    cell.cell_parameter_updated = true;
    basis.npwk_max = 9;
    kv.isk = {1, 0, 1};
    initialize(next_wavefunction);
    EXPECT_NEAR(radial_zero(), initial / 2.0, 1e-13);
    EXPECT_EQ(tools_created, 2);
    EXPECT_EQ(tools_destroyed, 1);
    EXPECT_EQ(tools_psi, &next_wavefunction);
    check_state(12, 9);
    // Exercise contraction as well as growth, keeping the extended radial range.
    cell.omega = 2.0;
    basis.npwk_max = 2;
    initialize(wavefunction);
    EXPECT_NEAR(radial_zero(), initial * 2.0, 1e-13);
    EXPECT_EQ(tools_created, 3);
    EXPECT_EQ(tools_destroyed, 2);
    check_state(12, 2);
}

TEST_F(OnsiteProjectorInitTest, UnchangedCellKeepsTableAndTools)
{
    const double initial = radial_zero();
    const std::complex<double>* buffer = projector.get_tab_atomic();
    initialize(wavefunction);
    EXPECT_DOUBLE_EQ(radial_zero(), initial);
    EXPECT_EQ(projector.get_tab_atomic(), buffer);
    EXPECT_EQ(tools_created, 1);
    EXPECT_EQ(tools_destroyed, 0);
    check_state(12, 4);
}

TEST_F(OnsiteProjectorInitTest, ConstantVolumeCellChangeStillRefreshesTools)
{
    const double initial = radial_zero();
    cell.cell_parameter_updated = true;
    basis.npwk_max = 7;
    initialize(wavefunction);
    EXPECT_DOUBLE_EQ(radial_zero(), initial);
    EXPECT_EQ(tools_created, 2);
    EXPECT_EQ(tools_destroyed, 1);
    check_state(12, 7);
}
}
