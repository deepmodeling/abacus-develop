#include "source_lcao/module_rt/snap_psb_half_tddft.h"

#include "source_base/ylm.h"
#include "source_cell/read_pp.h"
#include "source_cell/unitcell.h"
#include "source_io/module_hs/cal_r_overlap_r.h"
#include "../../lcao_nonlocal_info.h"

#ifdef __CUDA
#include "source_lcao/module_rt/kernels/snap_psibeta_gpu.h"
#include <cuda_runtime.h>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <vector>

SepPot::SepPot() = default;

SepPot::~SepPot()
{
    delete[] r;
    delete[] rv;
}

Sep_Cell::Sep_Cell() noexcept : ntype(0), omega(0.0), tpiba2(0.0)
{
}

Sep_Cell::~Sep_Cell() noexcept = default;

Magnetism::Magnetism()
{
    tot_mag = 0.0;
    abs_mag = 0.0;
}

Magnetism::~Magnetism()
{
}

UnitCell::UnitCell()
{
    itia2iat.create(1, 1);
}

UnitCell::~UnitCell()
{
    if (set_atom_flag)
    {
        delete[] atoms;
    }
}

void UnitCell::set_iat2iwt(const int& npol_in)
{
    npol = npol_in;
    iat2iwt.assign(nat, 0);
}

namespace
{
struct ComparisonStats
{
    double max_overlap_diff = 0.0;
    double max_position_diff = 0.0;
    double max_imag_abs = 0.0;
    double max_reference_abs = 0.0;
};

void print_comparison_stats(const char* label, const ComparisonStats& stats)
{
    std::cout << std::scientific << std::setprecision(12) << "psibeta " << label
              << ": max overlap error = " << stats.max_overlap_diff
              << ", max position error = " << stats.max_position_diff
              << ", max imaginary magnitude = " << stats.max_imag_abs
              << ", max reference magnitude = " << stats.max_reference_abs << std::defaultfloat << std::endl;
}

ComparisonStats compare_zero_vector_potential(const LCAO_Orbitals& orb,
                                              const UnitCell& ucell,
                                              cal_r_overlap_R& r_calculator,
                                              const int radial_grid_num,
                                              const int lebedev_grid_points)
{
    const ModuleBase::Vector3<double> R0(0.1, -0.2, 0.3);
    const ModuleBase::Vector3<double> R1(0.4, 0.2, -0.1);
    const ModuleBase::Vector3<double> zero_A(0.0, 0.0, 0.0);
    module_rt::SnapIntegrationOptions options;
    options.radial_grid_num = radial_grid_num;
    options.lebedev_grid_points = lebedev_grid_points;

    ComparisonStats stats;
    const auto* lcao_nl = dynamic_cast<const LCAONonlocalInfo*>(ucell.infoNL.get());
    EXPECT_NE(lcao_nl, nullptr);
    if (lcao_nl == nullptr)
    {
        return stats;
    }

    for (int L1 = 0; L1 <= orb.Phi[0].getLmax(); ++L1)
    {
        for (int N1 = 0; N1 < orb.Phi[0].getNchi(L1); ++N1)
        {
            for (int m1 = 0; m1 < 2 * L1 + 1; ++m1)
            {
                std::vector<std::vector<std::complex<double>>> grid_nlm;
                module_rt::snap_psibeta_half_tddft(orb,
                                                   lcao_nl->get_nonlocal(),
                                                   grid_nlm,
                                                   R1,
                                                   0,
                                                   L1,
                                                   m1,
                                                   N1,
                                                   R0,
                                                   0,
                                                   zero_A,
                                                   true,
                                                   options);

                std::vector<std::vector<double>> reference_nlm;
                r_calculator.get_psi_r_beta(ucell, reference_nlm, R1, 0, L1, m1, N1, R0, 0);

                EXPECT_EQ(grid_nlm.size(), 4);
                EXPECT_EQ(reference_nlm.size(), 4);
                if (grid_nlm.size() != 4 || reference_nlm.size() != 4)
                {
                    continue;
                }

                bool sizes_match = true;
                for (size_t dim = 0; dim < grid_nlm.size(); ++dim)
                {
                    EXPECT_EQ(grid_nlm[dim].size(), reference_nlm[dim].size());
                    sizes_match = sizes_match && (grid_nlm[dim].size() == reference_nlm[dim].size());
                }
                if (!sizes_match)
                {
                    continue;
                }

                for (size_t dim = 0; dim < grid_nlm.size(); ++dim)
                {
                    for (size_t i = 0; i < grid_nlm[dim].size(); ++i)
                    {
                        const double real_diff = std::abs(grid_nlm[dim][i].real() - reference_nlm[dim][i]);
                        if (dim == 0)
                        {
                            stats.max_overlap_diff = std::max(stats.max_overlap_diff, real_diff);
                        }
                        else
                        {
                            stats.max_position_diff = std::max(stats.max_position_diff, real_diff);
                        }
                        stats.max_imag_abs = std::max(stats.max_imag_abs, std::abs(grid_nlm[dim][i].imag()));
                        stats.max_reference_abs = std::max(stats.max_reference_abs, std::abs(reference_nlm[dim][i]));
                    }
                }
            }
        }
    }

    return stats;
}

#ifdef __CUDA
void compare_cpu_and_gpu_overlap(const LCAO_Orbitals& orb,
                                 const UnitCell& ucell,
                                 const Parallel_Orbitals& pv,
                                 const ModuleBase::Vector3<double>& vector_potential,
                                 const bool calc_r)
{
    int device_count = 0;
    const cudaError_t device_status = cudaGetDeviceCount(&device_count);
    if (device_status != cudaSuccess)
    {
        GTEST_SKIP() << "cudaGetDeviceCount failed: " << cudaGetErrorString(device_status);
    }
    if (device_count == 0)
    {
        GTEST_SKIP() << "cudaGetDeviceCount reported zero CUDA devices";
    }

    const ModuleBase::Vector3<double> R0(0.1, -0.2, 0.3);
    const ModuleBase::Vector3<double> R1(0.4, 0.2, -0.1);
    // The GPU production interface converts adjacent_tau from lattice units
    // to Cartesian coordinates, while the scalar CPU interface accepts R1
    // directly in Cartesian coordinates.
    const ModuleBase::Vector3<double> tau1(R1.x / ucell.lat0, R1.y / ucell.lat0, R1.z / ucell.lat0);
    const int nlm_dim = calc_r ? 4 : 1;

    AdjacentAtomInfo adjs;
    adjs.adj_num = 0;
    adjs.ntype.push_back(0);
    adjs.natom.push_back(0);
    adjs.adjacent_tau.push_back(tau1);
    adjs.box.push_back(ModuleBase::Vector3<int>(0, 0, 0));

    std::vector<std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>> gpu_nlm(
        1,
        std::vector<std::unordered_map<int, std::vector<std::complex<double>>>>(nlm_dim));

    const auto* lcao_nl = dynamic_cast<const LCAONonlocalInfo*>(ucell.infoNL.get());
    ASSERT_NE(lcao_nl, nullptr);
    module_rt::gpu::init_snap_psibeta_gpu();
    module_rt::gpu::snap_psibeta_atom_batch_gpu(orb,
                                                lcao_nl->get_nonlocal(),
                                                0,
                                                R0,
                                                vector_potential,
                                                adjs,
                                                &ucell,
                                                &pv,
                                                1,
                                                nlm_dim,
                                                gpu_nlm);

    const Atom& atom = ucell.atoms[0];
    double max_difference = 0.0;
    double max_cpu_abs = 0.0;
    for (int iw = 0; iw < atom.nw; ++iw)
    {
        std::vector<std::vector<std::complex<double>>> cpu_nlm;
        module_rt::snap_psibeta_half_tddft(orb,
                                           lcao_nl->get_nonlocal(),
                                           cpu_nlm,
                                           R1,
                                           0,
                                           atom.iw2l[iw],
                                           atom.iw2m[iw],
                                           atom.iw2n[iw],
                                           R0,
                                           0,
                                           vector_potential,
                                           calc_r);

        ASSERT_EQ(cpu_nlm.size(), static_cast<size_t>(nlm_dim));
        for (int dim = 0; dim < nlm_dim; ++dim)
        {
            const auto gpu_entry = gpu_nlm[0][dim].find(iw);
            ASSERT_NE(gpu_entry, gpu_nlm[0][dim].end());
            ASSERT_EQ(gpu_entry->second.size(), cpu_nlm[dim].size());
            for (size_t i = 0; i < cpu_nlm[dim].size(); ++i)
            {
                const double tolerance = 3.0e-14;
                const double difference = std::abs(gpu_entry->second[i] - cpu_nlm[dim][i]);
                max_difference = std::max(max_difference, difference);
                max_cpu_abs = std::max(max_cpu_abs, std::abs(cpu_nlm[dim][i]));
                EXPECT_LE(difference, tolerance)
                    << "iw = " << iw << ", dim = " << dim << ", projector component = " << i;
            }
        }
    }
    std::cout << std::scientific << std::setprecision(12)
              << "psibeta CPU/GPU calc_r=" << calc_r << ", A=(" << vector_potential.x << ", " << vector_potential.y
              << ", " << vector_potential.z << "): max difference = " << max_difference
              << ", max CPU magnitude = " << max_cpu_abs << std::defaultfloat << std::endl;
}
#endif

class SnapPsibetaHalfTddftTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ModuleBase::Ylm::set_coefficients();

        const std::string root = "../../../../../";
        const std::string orb_file = "tests/PP_ORB/Ti_gga_10au_100Ry_4s2p2d1f.orb";
        const std::string orbital_files[1] = {orb_file};

        std::ofstream ofs("snap_psibeta_half_tddft_test.log");
        // Keep the reciprocal transform spacing intentionally different from
        // the 0.01 Bohr real-space orbital grid.
        orb.init(ofs, 1, root, orbital_files, "", 3, 100.0, 0.013, 0.01, 30.0, false, 0, false, false, 0);

        ASSERT_EQ(orb.Phi[0].getLmax(), 3);
        ASSERT_EQ(orb.Phi[0].getNchi(0), 4);
        ASSERT_EQ(orb.Phi[0].getNchi(1), 2);
        ASSERT_EQ(orb.Phi[0].getNchi(2), 2);
        ASSERT_EQ(orb.Phi[0].getNchi(3), 1);

        build_ti_beta_projectors(root);
        initialize_r_overlap_reference();
    }

    void build_ti_beta_projectors(const std::string& root)
    {
        ucell.ntype = 1;
        ucell.nat = 1;
        ucell.lat0 = 2.0;
        ucell.atoms = new Atom[1];
        ucell.set_atom_flag = true;

        Atom& atom = ucell.atoms[0];
        atom.label = "Ti";
        atom.type = 0;
        atom.na = 1;
        atom.nwl = orb.Phi[0].getLmax();
        atom.l_nchi.resize(atom.nwl + 1);
        atom.nw = 0;
        for (int L = 0; L <= atom.nwl; ++L)
        {
            atom.l_nchi[L] = orb.Phi[0].getNchi(L);
            atom.nw += (2 * L + 1) * atom.l_nchi[L];
        }
        atom.tau.resize(1);
        atom.tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        atom.set_index();
        ucell.itia2iat(0, 0) = 0;
        ucell.set_iat2iwt(1);
        pv.set_serial(atom.nw, atom.nw);
        pv.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, atom.nw);

        Pseudopot_upf pseudo_reader;
        std::string pseudo_type = "auto";
        const int pseudo_error
            = pseudo_reader.init_pseudo_reader(root + "tests/PP_ORB/Ti_ONCV_PBE-1.0.upf", pseudo_type, atom.ncpp);
        ASSERT_EQ(pseudo_error, 0);
        ASSERT_EQ(pseudo_type, "upf201");
        ASSERT_EQ(atom.ncpp.psd, "Ti");
        ASSERT_EQ(atom.ncpp.pp_type, "NC");
        ASSERT_EQ(atom.ncpp.nbeta, 6);
        ASSERT_EQ(atom.ncpp.lll, std::vector<int>({0, 0, 1, 1, 2, 2}));
        const double pseudo_rcut = 15.0;
        pseudo_reader.complete_default(atom.ncpp, pseudo_rcut);
        ASSERT_EQ(atom.ncpp.nh, 18);
        ASSERT_EQ(atom.ncpp.jjj.size(), 6);

        auto* lcao_nl = new LCAONonlocalInfo();
        lcao_nl->get_nonlocal().nproj = new int[1];
        std::ofstream log("snap_psibeta_half_tddft_nonlocal.log");
        lcao_nl->get_nonlocal().Set_NonLocal(0,
                                            &atom,
                                            lcao_nl->get_nonlocal().nproj[0],
                                            orb.get_kmesh(),
                                            orb.get_dk(),
                                            orb.get_dr_uniform(),
                                            log,
                                            false,
                                            false,
                                            1);

        ASSERT_EQ(lcao_nl->get_nonlocal().nproj[0], 6);
        lcao_nl->get_nonlocal().nprojmax = lcao_nl->get_nonlocal().nproj[0];
        lcao_nl->get_nonlocal().rcutmax_Beta = lcao_nl->get_nonlocal().Beta[0].get_rcut_max();
        ucell.infoNL.reset(lcao_nl);
    }

    void initialize_r_overlap_reference()
    {
        r_calculator.init_nonlocal(ucell, pv, orb);
    }

    ComparisonStats compare_zero_vector_potential(const int radial_grid_num, const int lebedev_grid_points)
    {
        return ::compare_zero_vector_potential(orb, ucell, r_calculator, radial_grid_num, lebedev_grid_points);
    }

    LCAO_Orbitals orb;
    UnitCell ucell;
    Parallel_Orbitals pv;
    cal_r_overlap_R r_calculator;
};

class SnapPsibetaNonuniformHalfTddftTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        ModuleBase::Ylm::set_coefficients();

        const std::string root = "../../../../../";
        const std::string orbital_files[1] = {"tests/PP_ORB/Al_gga_10au_100Ry_3s3p2d.orb"};
        std::ofstream ofs("snap_psibeta_half_tddft_al_test.log");
        orb.init(ofs, 1, root, orbital_files, "", 2, 100.0, 0.017, 0.01, 30.0, false, 0, false, false, 0);

        ucell.ntype = 1;
        ucell.nat = 1;
        ucell.lat0 = 2.0;
        ucell.atoms = new Atom[1];
        ucell.set_atom_flag = true;

        Atom& atom = ucell.atoms[0];
        atom.label = "Al";
        atom.type = 0;
        atom.na = 1;
        atom.nwl = orb.Phi[0].getLmax();
        atom.l_nchi.resize(atom.nwl + 1);
        atom.nw = 0;
        for (int L = 0; L <= atom.nwl; ++L)
        {
            atom.l_nchi[L] = orb.Phi[0].getNchi(L);
            atom.nw += (2 * L + 1) * atom.l_nchi[L];
        }
        atom.tau.resize(1);
        atom.tau[0] = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
        atom.set_index();
        ucell.itia2iat(0, 0) = 0;
        ucell.set_iat2iwt(1);
        pv.set_serial(atom.nw, atom.nw);
        pv.set_atomic_trace(ucell.get_iat2iwt(), ucell.nat, atom.nw);

        Pseudopot_upf pseudo_reader;
        std::string pseudo_type = "auto";
        const int pseudo_error
            = pseudo_reader.init_pseudo_reader(root + "tests/PP_ORB/Al.pbe-rrkj.UPF", pseudo_type, atom.ncpp);
        ASSERT_EQ(pseudo_error, 0);
        ASSERT_EQ(atom.ncpp.psd, "Al");
        ASSERT_EQ(atom.ncpp.pp_type, "NC");
        ASSERT_EQ(atom.ncpp.nbeta, 4);
        pseudo_reader.complete_default(atom.ncpp, 15.0);

        auto* lcao_nl = new LCAONonlocalInfo();
        lcao_nl->get_nonlocal().nproj = new int[1];
        std::ofstream log("snap_psibeta_half_tddft_al_nonlocal.log");
        lcao_nl->get_nonlocal().Set_NonLocal(0,
                                            &atom,
                                            lcao_nl->get_nonlocal().nproj[0],
                                            orb.get_kmesh(),
                                            orb.get_dk(),
                                            orb.get_dr_uniform(),
                                            log,
                                            false,
                                            false,
                                            1);
        ASSERT_EQ(lcao_nl->get_nonlocal().nproj[0], 4);
        lcao_nl->get_nonlocal().nprojmax = lcao_nl->get_nonlocal().nproj[0];
        lcao_nl->get_nonlocal().rcutmax_Beta = lcao_nl->get_nonlocal().Beta[0].get_rcut_max();
        ucell.infoNL.reset(lcao_nl);

        r_calculator.init_nonlocal(ucell, pv, orb);
    }

    LCAO_Orbitals orb;
    UnitCell ucell;
    Parallel_Orbitals pv;
    cal_r_overlap_R r_calculator;
};
} // namespace

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialMatchesTwoCenterIntegral)
{
    const double overlap_tolerance = 5.0e-8;
    const double position_tolerance = 5.0e-8;
    const double imag_tolerance = 1.0e-14;
    const ComparisonStats stats = compare_zero_vector_potential(140, 110);
    print_comparison_stats("Ti 140x110", stats);

    EXPECT_LT(stats.max_overlap_diff, overlap_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_position_diff, position_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialDenseRadialGridMatchesTwoCenterIntegral)
{
    const double overlap_tolerance = 5.0e-8;
    const double position_tolerance = 6.0e-8;
    const double imag_tolerance = 1.0e-14;
    const ComparisonStats stats = compare_zero_vector_potential(280, 110);
    print_comparison_stats("Ti 280x110", stats);

    EXPECT_LT(stats.max_overlap_diff, overlap_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_position_diff, position_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaHalfTddftTest, ZeroVectorPotentialHighOrderGridMatchesTwoCenterIntegral)
{
    const double overlap_tolerance = 5.0e-8;
    const double position_tolerance = 5.0e-8;
    const double imag_tolerance = 1.0e-14;
    const ComparisonStats stats = compare_zero_vector_potential(140, 590);
    print_comparison_stats("Ti 140x590", stats);

    EXPECT_LT(stats.max_overlap_diff, overlap_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_position_diff, position_tolerance) << "max reference abs = " << stats.max_reference_abs;
    EXPECT_LT(stats.max_imag_abs, imag_tolerance) << "max reference abs = " << stats.max_reference_abs;
}

TEST_F(SnapPsibetaNonuniformHalfTddftTest, NonuniformAlProjectorMatchesTwoCenterIntegral)
{
    const auto* lcao_nl = dynamic_cast<const LCAONonlocalInfo*>(ucell.infoNL.get());
    ASSERT_NE(lcao_nl, nullptr);

    bool found_nonuniform_spacing = false;
    for (int ip = 0; ip < lcao_nl->get_nonlocal().nproj[0]; ++ip)
    {
        const auto& projector = lcao_nl->get_nonlocal().Beta[0].Proj[ip];
        ASSERT_GT(projector.getNr(), 2);
        const double first_spacing = projector.getRadial(1) - projector.getRadial(0);
        for (int ir = 2; ir < projector.getNr(); ++ir)
        {
            const double spacing = projector.getRadial(ir) - projector.getRadial(ir - 1);
            if (std::abs(spacing - first_spacing) > 1.0e-12)
            {
                found_nonuniform_spacing = true;
                break;
            }
        }
    }
    EXPECT_TRUE(found_nonuniform_spacing);

    const ComparisonStats default_grid = compare_zero_vector_potential(orb, ucell, r_calculator, 140, 110);
    const ComparisonStats dense_radial_grid = compare_zero_vector_potential(orb, ucell, r_calculator, 280, 110);
    const ComparisonStats dense_angular_grid = compare_zero_vector_potential(orb, ucell, r_calculator, 140, 590);
    print_comparison_stats("Al 140x110", default_grid);
    print_comparison_stats("Al 280x110", dense_radial_grid);
    print_comparison_stats("Al 140x590", dense_angular_grid);
    EXPECT_LT(default_grid.max_overlap_diff, 2.0e-3)
        << "max reference abs = " << default_grid.max_reference_abs;
    EXPECT_LT(default_grid.max_position_diff, 3.0e-3)
        << "max reference abs = " << default_grid.max_reference_abs;
    EXPECT_LT(dense_radial_grid.max_overlap_diff, 2.0e-5)
        << "max reference abs = " << dense_radial_grid.max_reference_abs;
    EXPECT_LT(dense_radial_grid.max_position_diff, 4.0e-5)
        << "max reference abs = " << dense_radial_grid.max_reference_abs;
    EXPECT_LT(dense_radial_grid.max_overlap_diff, 0.02 * default_grid.max_overlap_diff);
    EXPECT_LT(dense_radial_grid.max_position_diff, 0.02 * default_grid.max_position_diff);
    EXPECT_NEAR(dense_angular_grid.max_overlap_diff, default_grid.max_overlap_diff, 3.0e-9);
    EXPECT_NEAR(dense_angular_grid.max_position_diff, default_grid.max_position_diff, 1.0e-9);
    EXPECT_LT(default_grid.max_imag_abs, 1.0e-14)
        << "max reference abs = " << default_grid.max_reference_abs;
    EXPECT_LT(dense_radial_grid.max_imag_abs, 1.0e-14)
        << "max reference abs = " << dense_radial_grid.max_reference_abs;
    EXPECT_LT(dense_angular_grid.max_imag_abs, 1.0e-14)
        << "max reference abs = " << dense_angular_grid.max_reference_abs;
}

#ifdef __CUDA
TEST_F(SnapPsibetaHalfTddftTest, UniformTiCpuGpuOverlapWithoutPositionOperator)
{
    compare_cpu_and_gpu_overlap(orb, ucell, pv, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), false);
}

TEST_F(SnapPsibetaHalfTddftTest, UniformTiCpuGpuOverlapWithPositionOperator)
{
    compare_cpu_and_gpu_overlap(orb, ucell, pv, ModuleBase::Vector3<double>(0.03, -0.02, 0.01), true);
}

TEST_F(SnapPsibetaNonuniformHalfTddftTest, NonuniformAlCpuGpuOverlapWithoutPositionOperator)
{
    compare_cpu_and_gpu_overlap(orb, ucell, pv, ModuleBase::Vector3<double>(0.0, 0.0, 0.0), false);
}

TEST_F(SnapPsibetaNonuniformHalfTddftTest, NonuniformAlCpuGpuOverlapWithPositionOperator)
{
    compare_cpu_and_gpu_overlap(orb, ucell, pv, ModuleBase::Vector3<double>(0.03, -0.02, 0.01), true);
}
#endif
