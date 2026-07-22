#include "source_base/constants.h"
#include "source_base/math_ylmreal.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis_k.h"
#include "source_cell/unitcell.h"
#include "source_io/module_bessel/bessel_basis.h"
#include "source_io/module_bessel/numerical_basis.h"
#include "source_io/module_parameter/parameter.h"
#include "source_pw/module_pwdft/structure_factor.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace
{

using Complex = std::complex<double>;

constexpr double ecut_ry = 25.0;
constexpr double rcut_bohr = 2.0;
constexpr double tolerance = 1.0e-12;
constexpr int ik = 0;
const ModuleBase::Vector3<double> atom_tau(0.137, 0.211, 0.173);

TEST(SternheimerSIABPrimitives, DefaultPrimitiveEcutInheritsWavefunctionCutoff)
{
    Input_para& input = const_cast<Input_para&>(PARAM.inp);
    const std::string old_bessel_ecut = input.bessel_nao_ecut;
    const double old_ecutwfc = input.ecutwfc;
    const std::vector<double> old_rcuts = input.bessel_nao_rcuts;
    const int old_siab_lmax = input.sternheimer_siab_lmax;
    input.bessel_nao_ecut = "default";
    input.ecutwfc = 37.5;
    input.bessel_nao_rcuts = {8.0};
    input.sternheimer_siab_lmax = 2;

    const auto parameters = Numerical_Basis::siab_parameters_from_input(0, input.sternheimer_siab_lmax);
    EXPECT_DOUBLE_EQ(parameters.ecut_ry, 37.5);
    EXPECT_DOUBLE_EQ(parameters.rcut_bohr, 8.0);
    EXPECT_EQ(parameters.lmax, 2);

    input.bessel_nao_ecut = old_bessel_ecut;
    input.ecutwfc = old_ecutwfc;
    input.bessel_nao_rcuts = old_rcuts;
    input.sternheimer_siab_lmax = old_siab_lmax;
}

void expect_complex_near(const Complex& actual, const Complex& expected)
{
    const double scale = std::max(std::abs(actual), std::abs(expected));
    const double threshold = 1.0e-10 + 1.0e-10 * scale;
    EXPECT_NEAR(actual.real(), expected.real(), threshold);
    EXPECT_NEAR(actual.imag(), expected.imag(), threshold);
}

Complex inner_product(const std::vector<Complex>& left, const std::vector<Complex>& right, const double weight)
{
    EXPECT_EQ(left.size(), right.size());
    Complex local_result(0.0, 0.0);
    for (std::size_t i = 0; i < left.size(); ++i)
    {
        local_result += std::conj(left[i]) * right[i] * weight;
    }
#ifdef __MPI
    double local_values[2] = {local_result.real(), local_result.imag()};
    double global_values[2] = {0.0, 0.0};
    MPI_Allreduce(local_values, global_values, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return Complex(global_values[0], global_values[1]);
#else
    return local_result;
#endif
}

void fill_h_cell(UnitCell& ucell)
{
    ucell.lat0 = 8.0;
    ucell.latvec = ModuleBase::Matrix3(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    ucell.omega = 512.0;
    ucell.tpiba = ModuleBase::TWO_PI / ucell.lat0;
    ucell.ntype = 1;
    ucell.nat = 1;
    ucell.lmax = 1;
    ucell.lmaxmax = 1;
    ucell.nmax = 1;
    ucell.atoms = new Atom[1];
    ucell.atoms[0].type = 0;
    ucell.atoms[0].label = "H";
    ucell.atoms[0].na = 1;
    ucell.atoms[0].nwl = 1;
    ucell.atoms[0].l_nchi = {1, 1};
    ucell.atoms[0].tau = {atom_tau};
}

std::vector<Complex> reciprocal_primitive(const Bessel_Basis& bessel,
                                          const ModulePW::PW_Basis_K& pw,
                                          const UnitCell& ucell,
                                          const int l,
                                          const int conventional_m,
                                          const int ie)
{
    const int npw = pw.npwk[ik];
    std::vector<ModuleBase::Vector3<double>> gk(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        gk[ig] = pw.getgpluskcar(ik, ig) * ucell.tpiba;
    }

    const int target_lmax = std::max(ucell.lmax, l);
    const int total_lm = (target_lmax + 1) * (target_lmax + 1);
    ModuleBase::matrix ylm(total_lm, npw);
    ModuleBase::YlmReal::Ylm_Real(total_lm, npw, gk.data(), ylm);

    const int abacus_m = Numerical_Basis::siab_abacus_m(conventional_m);
    const int lm = l * l + abacus_m;
    const Complex lphase = (4.0 * ModuleBase::PI / std::sqrt(ucell.omega))
                           * std::pow(ModuleBase::IMAG_UNIT, -l);

    std::vector<Complex> values(npw);
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> gdirect = pw.getgdirect(ik, ig);
        const Complex sk = std::exp(Complex(0.0, -ModuleBase::TWO_PI * (gdirect * atom_tau)));
        values[ig] = lphase * sk * ylm(lm, ig)
                     * bessel.Polynomial_Interpolation2(l, ie, gk[ig].norm());
    }
    return values;
}

std::vector<Complex> physical_complex_grid(const ModulePW::PW_Basis_K& pw,
                                           const UnitCell& ucell,
                                           const std::vector<Complex>& reciprocal)
{
    EXPECT_FALSE(pw.gamma_only);
    std::vector<Complex> grid(std::max(1, pw.nrxx), Complex(0.0, 0.0));
    pw.recip2real(reciprocal.data(), grid.data(), ik);
    grid.resize(pw.nrxx);
    const double scale = 1.0 / std::sqrt(ucell.omega);
    for (Complex& value : grid)
    {
        value *= scale;
    }
    return grid;
}

std::vector<Complex> deterministic_complex_pw(const ModulePW::PW_Basis_K& pw)
{
    std::vector<Complex> values(pw.npwk[ik]);
    for (int ig = 0; ig < pw.npwk[ik]; ++ig)
    {
        const ModuleBase::Vector3<double> g = pw.getgdirect(ik, ig);
        const double denominator = 1.0 + g.norm2();
        values[ig] = Complex((1.0 + 0.17 * g.x - 0.11 * g.z) / denominator,
                             (0.23 + 0.13 * g.y + 0.07 * g.z) / denominator);
    }
    return values;
}

void initialize_reference_bessel(Bessel_Basis& bessel,
                                 const UnitCell& ucell,
                                 const Numerical_Basis::SIABPrimitiveParameters& parameters)
{
    const int target_lmax = parameters.lmax < 0 ? ucell.lmax : parameters.lmax;
    bessel.init(false,
                parameters.ecut_ry,
                ucell.ntype,
                target_lmax,
                parameters.smooth,
                parameters.sigma,
                parameters.rcut_bohr,
                parameters.tolerance,
                ucell);
}

class SternheimerSIABPrimitivesTest : public ::testing::Test
{
  protected:
    UnitCell ucell;
    ModulePW::PW_Basis_K gamma_pw;
    ModulePW::PW_Basis_K full_pw;
    Structure_Factor structure_factor;

    void SetUp() override
    {
        fill_h_cell(ucell);
        const ModuleBase::Vector3<double> gamma(0.0, 0.0, 0.0);
#ifdef __MPI
        int rank = 0;
        int size = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        gamma_pw.initmpi(size, rank, MPI_COMM_WORLD);
        full_pw.initmpi(size, rank, MPI_COMM_WORLD);
#endif
        gamma_pw.initgrids(ucell.lat0, ucell.latvec, 4.0 * ecut_ry);
        gamma_pw.initparameters(true, ecut_ry, 1, &gamma, 2, true);
        gamma_pw.setuptransform();
        gamma_pw.collect_local_pw();

        full_pw.initgrids(ucell.lat0, ucell.latvec, gamma_pw.nx, gamma_pw.ny, gamma_pw.nz);
        full_pw.initparameters(false, ecut_ry, 1, &gamma, 2, true);
        full_pw.setuptransform();
        full_pw.collect_local_pw();

        ASSERT_TRUE(gamma_pw.gamma_only);
        ASSERT_FALSE(full_pw.gamma_only);
        ASSERT_EQ(gamma_pw.nx, full_pw.nx);
        ASSERT_EQ(gamma_pw.ny, full_pw.ny);
        ASSERT_EQ(gamma_pw.nz, full_pw.nz);
        ASSERT_EQ(gamma_pw.nxyz, full_pw.nxyz);
        ASSERT_EQ(gamma_pw.nplane, full_pw.nplane);
        ASSERT_EQ(gamma_pw.startz_current, full_pw.startz_current);
        ASSERT_EQ(gamma_pw.nrxx, full_pw.nrxx);
    }

    void TearDown() override
    {
        delete[] ucell.atoms;
        ucell.atoms = nullptr;
    }
};

} // namespace

TEST(SternheimerSIABPrimitiveMapping, ConvertsConventionalMToABACUSOrdering)
{
    EXPECT_EQ(Numerical_Basis::siab_abacus_m(0), 0);
    EXPECT_EQ(Numerical_Basis::siab_abacus_m(1), 1);
    EXPECT_EQ(Numerical_Basis::siab_abacus_m(-1), 2);
    EXPECT_EQ(Numerical_Basis::siab_abacus_m(2), 3);
    EXPECT_EQ(Numerical_Basis::siab_abacus_m(-2), 4);
}

TEST_F(SternheimerSIABPrimitivesTest, OrdersBlocksAndPreservesFullComplexPWInnerProductsOnGrid)
{
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
    };
    Numerical_Basis numerical_basis;
    const auto blocks
        = numerical_basis.siab_primitive_grid_values(ik, &full_pw, structure_factor, ucell, parameters);

    const int nprimitive = static_cast<int>(std::sqrt(ecut_ry) * rcut_bohr / ModuleBase::PI);
    ASSERT_EQ(nprimitive, 3);
    ASSERT_EQ(blocks.size(), 4);

    EXPECT_EQ(blocks[0].type_index, 0);
    EXPECT_EQ(blocks[0].element, "H");
    EXPECT_EQ(blocks[0].atom_index, 0);
    EXPECT_EQ(blocks[0].l, 0);
    EXPECT_EQ(blocks[0].m, 0);
    EXPECT_EQ(blocks[0].n_primitive, nprimitive);
    EXPECT_EQ(blocks[0].offset, 0);

    EXPECT_EQ(blocks[1].l, 1);
    EXPECT_EQ(blocks[1].m, -1);
    EXPECT_EQ(blocks[1].offset, nprimitive);
    EXPECT_EQ(blocks[2].l, 1);
    EXPECT_EQ(blocks[2].m, 0);
    EXPECT_EQ(blocks[2].offset, 2 * nprimitive);
    EXPECT_EQ(blocks[3].l, 1);
    EXPECT_EQ(blocks[3].m, 1);
    EXPECT_EQ(blocks[3].offset, 3 * nprimitive);

    Bessel_Basis bessel;
    initialize_reference_bessel(bessel, ucell, parameters);

    const double delta_omega = ucell.omega / static_cast<double>(full_pw.nxyz);
    for (const auto& block : blocks)
    {
        ASSERT_EQ(block.n_primitive, nprimitive);
        ASSERT_EQ(block.values.size(), static_cast<std::size_t>(nprimitive));
        for (const auto& values : block.values)
        {
            EXPECT_EQ(values.size(), static_cast<std::size_t>(full_pw.nrxx));
        }

        for (int ie_left = 0; ie_left < nprimitive; ++ie_left)
        {
            const auto reciprocal_left
                = reciprocal_primitive(bessel, full_pw, ucell, block.l, block.m, ie_left);
            for (int ie_right = 0; ie_right < nprimitive; ++ie_right)
            {
                const auto reciprocal_right
                    = reciprocal_primitive(bessel, full_pw, ucell, block.l, block.m, ie_right);
                const Complex reciprocal_overlap = inner_product(reciprocal_left, reciprocal_right, 1.0);
                const Complex grid_overlap
                    = inner_product(block.values[ie_left], block.values[ie_right], delta_omega);
                expect_complex_near(grid_overlap, reciprocal_overlap);
            }
        }
    }
}

TEST_F(SternheimerSIABPrimitivesTest, ExposesTheSamePhysicallyNormalizedReciprocalPrimitives)
{
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
        1,
    };
    Numerical_Basis numerical_basis;
    const auto blocks
        = numerical_basis.siab_primitive_reciprocal_values(ik, &full_pw, structure_factor, ucell, parameters);
    Bessel_Basis bessel;
    initialize_reference_bessel(bessel, ucell, parameters);
    const int nprimitive = static_cast<int>(std::sqrt(ecut_ry) * rcut_bohr / ModuleBase::PI);
    ASSERT_EQ(blocks.size(), 4);
    for (const auto& block: blocks)
    {
        ASSERT_EQ(block.values.size(), static_cast<std::size_t>(nprimitive));
        for (int ie = 0; ie != nprimitive; ++ie)
        {
            const auto expected = reciprocal_primitive(bessel,
                                                       full_pw,
                                                       ucell,
                                                       block.l,
                                                       block.m,
                                                       ie);
            ASSERT_EQ(block.values[static_cast<std::size_t>(ie)].size(), expected.size());
            for (std::size_t ig = 0; ig != expected.size(); ++ig)
            {
                expect_complex_near(block.values[static_cast<std::size_t>(ie)][ig], expected[ig]);
            }
        }
    }
}

TEST_F(SternheimerSIABPrimitivesTest, RecoversNormalizedCoefficientsFromPhysicalGridValues)
{
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
        1,
    };
    Bessel_Basis bessel;
    initialize_reference_bessel(bessel, ucell, parameters);
    const auto expected = reciprocal_primitive(bessel, full_pw, ucell, 1, -1, 0);
    const auto physical_grid = physical_complex_grid(full_pw, ucell, expected);
    std::vector<Complex> recovered(static_cast<std::size_t>(full_pw.npwk[ik]), Complex(0.0, 0.0));
    full_pw.real2recip(physical_grid.data(), recovered.data(), ik);
    for (Complex& value: recovered)
    {
        value *= std::sqrt(ucell.omega);
    }
    ASSERT_EQ(recovered.size(), expected.size());
    for (std::size_t ig = 0; ig != expected.size(); ++ig)
    {
        expect_complex_near(recovered[ig], expected[ig]);
    }
}

TEST_F(SternheimerSIABPrimitivesTest, UsesGlobalAtomIndicesAcrossSpecies)
{
    delete[] ucell.atoms;
    ucell.ntype = 2;
    ucell.nat = 2;
    ucell.lmax = 0;
    ucell.lmaxmax = 0;
    ucell.nmax = 1;
    ucell.atoms = new Atom[2];
    for (int type = 0; type != 2; ++type)
    {
        ucell.atoms[type].type = type;
        ucell.atoms[type].label = type == 0 ? "H" : "H_empty";
        ucell.atoms[type].na = 1;
        ucell.atoms[type].nwl = 0;
        ucell.atoms[type].l_nchi = {1};
        ucell.atoms[type].tau = {atom_tau};
    }
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
        0,
    };
    Numerical_Basis numerical_basis;
    const auto blocks
        = numerical_basis.siab_primitive_reciprocal_values(ik, &full_pw, structure_factor, ucell, parameters);
    ASSERT_EQ(blocks.size(), 2);
    EXPECT_EQ(blocks[0].element, "H");
    EXPECT_EQ(blocks[0].atom_index, 0);
    EXPECT_EQ(blocks[1].element, "H_empty");
    EXPECT_EQ(blocks[1].atom_index, 1);
}

TEST_F(SternheimerSIABPrimitivesTest, ExplicitLmaxExtendsPrimitivesBeyondInputOrbitalChannels)
{
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
        2,
    };
    ASSERT_EQ(ucell.atoms[0].nwl, 1);

    Numerical_Basis numerical_basis;
    const auto blocks
        = numerical_basis.siab_primitive_grid_values(ik, &full_pw, structure_factor, ucell, parameters);

    const int nprimitive = static_cast<int>(std::sqrt(ecut_ry) * rcut_bohr / ModuleBase::PI);
    ASSERT_EQ(blocks.size(), 9);
    for (int magnetic_offset = 0; magnetic_offset != 5; ++magnetic_offset)
    {
        const auto& block = blocks[static_cast<std::size_t>(4 + magnetic_offset)];
        EXPECT_EQ(block.l, 2);
        EXPECT_EQ(block.m, magnetic_offset - 2);
        EXPECT_EQ(block.n_primitive, nprimitive);
        EXPECT_EQ(block.offset, (4 + magnetic_offset) * nprimitive);
    }

    Bessel_Basis bessel;
    initialize_reference_bessel(bessel, ucell, parameters);
    const auto& d_block = blocks.back();
    const auto reciprocal_reference
        = reciprocal_primitive(bessel, full_pw, ucell, d_block.l, d_block.m, 0);
    const auto grid_reference = physical_complex_grid(full_pw, ucell, reciprocal_reference);
    ASSERT_EQ(d_block.values[0].size(), grid_reference.size());
    for (std::size_t ir = 0; ir != grid_reference.size(); ++ir)
    {
        expect_complex_near(d_block.values[0][ir], grid_reference[ir]);
    }
}

TEST_F(SternheimerSIABPrimitivesTest, GammaCompressedOddLMatchesIndependentFullComplexGrid)
{
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
    };
    Numerical_Basis numerical_basis;
    const auto gamma_blocks
        = numerical_basis.siab_primitive_grid_values(ik, &gamma_pw, structure_factor, ucell, parameters);

    ASSERT_EQ(gamma_blocks.size(), 4);
    const auto& odd_l_block = gamma_blocks[3];
    ASSERT_EQ(odd_l_block.l, 1);
    ASSERT_EQ(odd_l_block.m, 1);
    ASSERT_FALSE(odd_l_block.values.empty());

    Bessel_Basis bessel;
    initialize_reference_bessel(bessel, ucell, parameters);
    const auto reciprocal_reference
        = reciprocal_primitive(bessel, full_pw, ucell, odd_l_block.l, odd_l_block.m, 0);
    const auto grid_reference = physical_complex_grid(full_pw, ucell, reciprocal_reference);

    ASSERT_EQ(odd_l_block.values[0].size(), grid_reference.size());
    for (std::size_t ir = 0; ir < grid_reference.size(); ++ir)
    {
        expect_complex_near(odd_l_block.values[0][ir], grid_reference[ir]);
    }
}

TEST_F(SternheimerSIABPrimitivesTest, FullComplexReciprocalBpsiMatchesGridOverlap)
{
    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
    };
    Numerical_Basis numerical_basis;
    const auto blocks
        = numerical_basis.siab_primitive_grid_values(ik, &full_pw, structure_factor, ucell, parameters);

    ASSERT_EQ(blocks.size(), 4);
    const auto& odd_l_block = blocks[3];
    ASSERT_EQ(odd_l_block.l, 1);
    ASSERT_EQ(odd_l_block.m, 1);

    Bessel_Basis bessel;
    initialize_reference_bessel(bessel, ucell, parameters);
    const auto reciprocal_b
        = reciprocal_primitive(bessel, full_pw, ucell, odd_l_block.l, odd_l_block.m, 0);
    const auto reciprocal_psi = deterministic_complex_pw(full_pw);
    const auto grid_psi = physical_complex_grid(full_pw, ucell, reciprocal_psi);

    const Complex reciprocal_overlap = inner_product(reciprocal_b, reciprocal_psi, 1.0);
    const double delta_omega = ucell.omega / static_cast<double>(full_pw.nxyz);
    const Complex grid_overlap = inner_product(odd_l_block.values[0], grid_psi, delta_omega);
    expect_complex_near(grid_overlap, reciprocal_overlap);
}

TEST_F(SternheimerSIABPrimitivesTest, ReinitializesWhenParametersOrCellBasisShapeChange)
{
    const Numerical_Basis::SIABPrimitiveParameters first_parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
    };
    const Numerical_Basis::SIABPrimitiveParameters second_parameters{
        36.0,
        3.0,
        false,
        0.1,
        tolerance,
    };

    Numerical_Basis numerical_basis;
    const auto first_blocks
        = numerical_basis.siab_primitive_grid_values(ik, &full_pw, structure_factor, ucell, first_parameters);
    ASSERT_EQ(first_blocks.size(), 4);
    EXPECT_EQ(first_blocks[0].n_primitive, 3);

    const auto second_blocks
        = numerical_basis.siab_primitive_grid_values(ik, &full_pw, structure_factor, ucell, second_parameters);
    ASSERT_EQ(second_blocks.size(), 4);
    EXPECT_EQ(second_blocks[0].n_primitive, 5);
    EXPECT_EQ(second_blocks[3].offset, 15);

    ucell.lmax = 2;
    ucell.lmaxmax = 2;
    ucell.atoms[0].nwl = 2;
    ucell.atoms[0].l_nchi = {1, 1, 1};
    const auto reshaped_blocks
        = numerical_basis.siab_primitive_grid_values(ik, &full_pw, structure_factor, ucell, second_parameters);
    ASSERT_EQ(reshaped_blocks.size(), 9);
    EXPECT_EQ(reshaped_blocks[4].l, 2);
    EXPECT_EQ(reshaped_blocks[4].m, -2);
    EXPECT_EQ(reshaped_blocks[8].l, 2);
    EXPECT_EQ(reshaped_blocks[8].m, 2);
    EXPECT_EQ(reshaped_blocks[8].n_primitive, 5);
    EXPECT_EQ(reshaped_blocks[8].offset, 40);
}

TEST(SternheimerSIABPrimitivesMPI, EmptyLocalPWStillParticipatesInGridTransform)
{
#ifndef __MPI
    GTEST_SKIP();
#else
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2)
    {
        GTEST_SKIP();
    }

    UnitCell ucell;
    fill_h_cell(ucell);
    const ModuleBase::Vector3<double> gamma(0.0, 0.0, 0.0);
    ModulePW::PW_Basis_K tiny_pw;
    tiny_pw.initmpi(size, rank, MPI_COMM_WORLD);
    tiny_pw.initgrids(ucell.lat0, ucell.latvec, 4.0 * ecut_ry);
    tiny_pw.initparameters(false, 4.0, 1, &gamma, 2, true);
    tiny_pw.setuptransform();
    tiny_pw.collect_local_pw();

    ASSERT_GT(tiny_pw.npwk[ik], 0);
    if (rank == size - 1)
    {
        // setupIndGk rejects empty ranks, so retain its FFT plan and emulate a tighter local cutoff.
        tiny_pw.npwk[ik] = 0;
    }
    const int local_empty = tiny_pw.npwk[ik] == 0 ? 1 : 0;
    int empty_ranks = 0;
    MPI_Allreduce(&local_empty, &empty_ranks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ASSERT_EQ(empty_ranks, 1);
    if (local_empty)
    {
        ASSERT_GT(tiny_pw.nrxx, 0);
    }

    const Numerical_Basis::SIABPrimitiveParameters parameters{
        ecut_ry,
        rcut_bohr,
        false,
        0.1,
        tolerance,
    };
    Structure_Factor structure_factor;
    Numerical_Basis numerical_basis;
    const auto blocks
        = numerical_basis.siab_primitive_grid_values(ik, &tiny_pw, structure_factor, ucell, parameters);

    ASSERT_EQ(blocks.size(), 4);
    for (const auto& block : blocks)
    {
        ASSERT_EQ(block.values.size(), 3);
        for (const auto& values : block.values)
        {
            EXPECT_EQ(values.size(), static_cast<std::size_t>(tiny_pw.nrxx));
        }
    }
    if (local_empty)
    {
        double local_grid_norm = 0.0;
        for (const Complex& value : blocks[0].values[0])
        {
            local_grid_norm += std::norm(value);
        }
        EXPECT_GT(local_grid_norm, 0.0);
    }

    delete[] ucell.atoms;
    ucell.atoms = nullptr;
#endif
}

Structure_Factor::Structure_Factor() = default;
Structure_Factor::~Structure_Factor() = default;

std::complex<double>* Structure_Factor::get_sk(const int ik_in,
                                               const int it,
                                               const int ia,
                                               const ModulePW::PW_Basis_K* wfc_basis) const
{
    EXPECT_GE(it, 0);
    EXPECT_GE(ia, 0);
    const int npw = wfc_basis->npwk[ik_in];
    std::complex<double>* sk = new std::complex<double>[npw];
    for (int ig = 0; ig < npw; ++ig)
    {
        const ModuleBase::Vector3<double> gdirect = wfc_basis->getgdirect(ik_in, ig);
        const ModuleBase::Vector3<double> shifted_tau
            = atom_tau + ModuleBase::Vector3<double>(0.013 * it, 0.017 * ia, 0.0);
        sk[ig] = std::exp(Complex(0.0, -ModuleBase::TWO_PI * (gdirect * shifted_tau)));
    }
    return sk;
}

Atom::Atom() = default;
Atom::~Atom() = default;
Atom_pseudo::Atom_pseudo() = default;
Atom_pseudo::~Atom_pseudo() = default;
pseudo::pseudo() = default;
pseudo::~pseudo() = default;
UnitCell::UnitCell() = default;
UnitCell::~UnitCell() = default;
Magnetism::Magnetism() = default;
Magnetism::~Magnetism() = default;
SepPot::SepPot() = default;
SepPot::~SepPot() = default;
Sep_Cell::Sep_Cell() noexcept = default;
Sep_Cell::~Sep_Cell() noexcept = default;

#ifdef __MPI
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
#endif
