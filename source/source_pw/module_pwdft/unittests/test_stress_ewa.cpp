#include "source_pw/module_pwdft/stress_func.h"
#include "source_base/constants.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include "source_hamilt/module_ewald/h_ewald_pw.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

// Parameter explicitly grants this test helper access; no access-control macros.
class TestParameters
{
  public:
    static void gamma_only(const bool enabled)
    {
        Parameter& parameters = PARAM;
        parameters.sys.gamma_only_pw = enabled;
    }
};

namespace
{
int process_count = 1;
int process_rank = 0;

double to_kbar()
{
    return ModuleBase::RYDBERG_SI / std::pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
}

double choose_alpha(const UnitCell& cell, const ModulePW::PW_Basis& basis)
{
    double charge = 0.0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        charge += cell.atoms[it].na * cell.atoms[it].ncpp.zv;
    }
    double alpha = 2.9;
    double upper_bound = 0.0;
    do
    {
        alpha -= 0.1;
        const double prefactor_argument = 2 * alpha / ModuleBase::TWO_PI;
        const double prefactor = std::sqrt(prefactor_argument);
        const double cutoff_argument = cell.tpiba2 * basis.ggecut / (4 * alpha);
        const double cutoff_root = std::sqrt(cutoff_argument);
        upper_bound = ModuleBase::e2 * charge * charge * prefactor * std::erfc(cutoff_root);
    } while (upper_bound > 1.0e-7);
    return alpha;
}

// Independent serial reference: the former per-atom reciprocal formula and
// the real-space pair virial, before any crystal symmetry is applied.
ModuleBase::matrix reference_stress(const UnitCell& cell,
                                    const ModulePW::PW_Basis& basis,
                                    const double reciprocal_weight)
{
    ModuleBase::matrix sigma(3, 3);
    const double alpha = choose_alpha(cell, basis);
    double charge = 0.0;
    for (int it = 0; it < cell.ntype; ++it)
    {
        charge += cell.atoms[it].na * cell.atoms[it].ncpp.zv;
    }
    if (basis.ig_gge0 >= 0)
    {
        const double charge_density = charge / cell.omega;
        const double background = ModuleBase::TWO_PI * ModuleBase::e2 / (4 * alpha)
                                  * std::pow(charge_density, 2);
        for (int i = 0; i < 3; ++i)
        {
            sigma(i, i) -= background;
        }
    }
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        if (ig == basis.ig_gge0)
        {
            continue;
        }
        std::complex<double> q(0.0, 0.0);
        for (int it = 0; it < cell.ntype; ++it)
        {
            for (const auto& tau : cell.atoms[it].tau)
            {
                const double phase = ModuleBase::TWO_PI * (basis.gcar[ig] * tau);
                const double real = std::cos(phase);
                const double imag = std::sin(phase);
                q += cell.atoms[it].ncpp.zv * std::complex<double>(real, imag);
            }
        }
        const double g2 = basis.gg[ig] * cell.tpiba2;
        const double a = g2 / (4 * alpha);
        const double exponent = -a;
        const double damping = std::exp(exponent);
        const std::complex<double> charge_density = q / cell.omega;
        const double energy_density = reciprocal_weight * ModuleBase::TWO_PI * ModuleBase::e2
                                      * damping / g2 * std::norm(charge_density);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                sigma(i, j) += energy_density
                               * ((i == j ? 1.0 : 0.0)
                                  - 2 * cell.tpiba2 * basis.gcar[ig][i] * basis.gcar[ig][j] / g2 * (a + 1));
            }
        }
    }
    if (basis.ig_gge0 >= 0)
    {
        const double sqrt_alpha = std::sqrt(alpha);
        const double real_prefactor_argument = 8 * alpha / ModuleBase::TWO_PI;
        const double real_prefactor = std::sqrt(real_prefactor_argument);
        const double rmax = 4 / sqrt_alpha / cell.lat0;
        const int mxr = H_Ewald_pw::estimate_mxr(rmax, cell.G);
        std::vector<int> order(mxr);
        std::vector<double> distance2(mxr);
        std::vector<ModuleBase::Vector3<double>> images(mxr);
        for (int a = 0; a < cell.nat; ++a)
        {
            for (int b = 0; b < cell.nat; ++b)
            {
                const int it = cell.iat2it[a];
                const int jt = cell.iat2it[b];
                const auto displacement = cell.atoms[it].tau[cell.iat2ia[a]] - cell.atoms[jt].tau[cell.iat2ia[b]];
                int count = 0;
                H_Ewald_pw::rgen(displacement, rmax, order.data(), cell.latvec, cell.G,
                                 images.data(), distance2.data(), mxr, count, 0);
                for (int n = 0; n < count; ++n)
                {
                    const double r = std::sqrt(distance2[n]) * cell.lat0;
                    const double erfc_argument = sqrt_alpha * r;
                    const double exponent = -alpha * r * r;
                    const double coefficient = ModuleBase::e2 / (2 * cell.omega)
                                               * cell.atoms[it].ncpp.zv * cell.atoms[jt].ncpp.zv
                                               * (std::erfc(erfc_argument)
                                                  + r * real_prefactor * std::exp(exponent))
                                               / (r * r * r) * cell.lat0 * cell.lat0;
                    for (int i = 0; i < 3; ++i)
                    {
                        for (int j = 0; j < 3; ++j)
                        {
                            sigma(i, j) += coefficient * images[n][i] * images[n][j];
                        }
                    }
                }
            }
        }
    }
    Parallel_Reduce::reduce_pool(sigma.c, 9);
    return sigma;
}

void expect_tensor(const ModuleBase::matrix& actual, const ModuleBase::matrix& expected, const double tolerance)
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            const double actual_kbar = actual(i, j) * to_kbar();
            const double expected_kbar = expected(i, j) * to_kbar();
            EXPECT_NEAR(actual_kbar, expected_kbar, tolerance) << i << "," << j;
            EXPECT_DOUBLE_EQ(actual(i, j), actual(j, i));
        }
    }
}

ModuleBase::Matrix3 strained_lattice(const ModuleBase::Matrix3& lattice,
                                      const int i,
                                      const int j,
                                      const double step)
{
    double strain[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const double increment = i == j ? step : step / 2;
    strain[i][j] += increment;
    if (i != j)
    {
        strain[j][i] += increment;
    }
    const ModuleBase::Matrix3 deformation(strain[0][0], strain[0][1], strain[0][2],
                                          strain[1][0], strain[1][1], strain[1][2],
                                          strain[2][0], strain[2][1], strain[2][2]);
    const ModuleBase::Matrix3 transpose = deformation.Transpose();
    return lattice * transpose;
}

class StressEwaTest : public testing::Test
{
  protected:
    UnitCell cell;
    std::vector<Atom> atoms;
    ModulePW::PW_Basis basis;
    Structure_Factor sf;
    Parallel_Grid grid;
    Stress_Func<double> stress;

    void SetUp() override
    {
        TestParameters::gamma_only(false);
        cell.lat0 = 8.0;
        cell.latvec = ModuleBase::Matrix3(1.0, 0.13, 0.04, 0.08, 1.1, 0.17, 0.03, 0.07, 0.95);
        update_geometry(cell.latvec);
        atoms.resize(2);
        cell.atoms = atoms.data();
        cell.ntype = 2;
        cell.nat = 3;
        atoms[0].na = 2;
        atoms[0].ncpp.zv = 1.2;
        atoms[0].taud = {{0.13, 0.24, 0.32}, {0.47, 0.56, 0.68}};
        atoms[1].na = 1;
        atoms[1].ncpp.zv = 3.7;
        atoms[1].taud = {{0.79, 0.18, 0.43}};
        update_positions();
        cell.set_iat2itia();
        basis.initmpi(process_count, process_rank, POOL_WORLD);
        basis.initgrids(cell.lat0, cell.latvec, 36, 36, 36);
        basis.initparameters(false, 80.0);
        basis.setuptransform();
        basis.collect_local_pw();
        sf.set(&basis, 0);
        sf.setup(&cell, grid, &basis, false);
    }

    void TearDown() override
    {
        TestParameters::gamma_only(false);
        // UnitCell owns its index maps; the vector owns the mock atoms.
    }

    void update_positions()
    {
        for (auto& atom : atoms)
        {
            atom.tau.resize(atom.na);
            for (int i = 0; i < atom.na; ++i)
            {
                atom.tau[i] = atom.taud[i] * cell.latvec;
            }
        }
    }

    void update_geometry(const ModuleBase::Matrix3& lattice)
    {
        cell.latvec = lattice;
        cell.GT = lattice.Inverse();
        cell.G = cell.GT.Transpose();
        const double determinant = lattice.Det();
        cell.omega = std::abs(determinant) * std::pow(cell.lat0, 3);
        cell.tpiba = ModuleBase::TWO_PI / cell.lat0;
        cell.tpiba2 = cell.tpiba * cell.tpiba;
        update_positions();
        for (int ig = 0; ig < basis.npw; ++ig)
        {
            basis.gcar[ig] = basis.gdirect[ig] * cell.G;
            basis.gg[ig] = basis.gcar[ig].norm2();
        }
    }

    ModuleBase::matrix calculate(const bool is_pw)
    {
        ModuleBase::matrix sigma(3, 3);
        stress.stress_ewa(cell, sigma, &basis, sf.strucFac, is_pw);
        return sigma;
    }
};

TEST_F(StressEwaTest, DirectSumChargesNormalizationAndFullTensor)
{
    const auto expected = reference_stress(cell, basis, 1.0);
    const auto actual = calculate(false);
    expect_tensor(actual, expected, 1.0e-8);
    const double shear_kbar = expected(0, 1) * to_kbar();
    EXPECT_GT(std::abs(shear_kbar), 1.0e-3);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        for (int it = 0; it < cell.ntype; ++it)
        {
            std::complex<double> negative_g(0, 0);
            for (const auto& tau : atoms[it].tau)
            {
                const double phase = ModuleBase::TWO_PI * (basis.gcar[ig] * tau);
                const double real = std::cos(phase);
                const double imag = std::sin(phase);
                negative_g += std::complex<double>(real, imag);
            }
            const std::complex<double> sf_value = sf.strucFac(it, ig);
            const std::complex<double> difference = std::conj(sf_value) - negative_g;
            EXPECT_NEAR(std::abs(difference), 0.0, 1.0e-12);
        }
        if (ig == basis.ig_gge0)
        {
            const std::complex<double> first_error = sf.strucFac(0, ig) - 2.0;
            const std::complex<double> second_error = sf.strucFac(1, ig) - 1.0;
            EXPECT_NEAR(std::abs(first_error), 0.0, 1.0e-14);
            EXPECT_NEAR(std::abs(second_error), 0.0, 1.0e-14);
        }
    }
    // A supplied SF must be consumed. In particular, S(0) must be skipped.
    if (basis.ig_gge0 >= 0)
    {
        sf.strucFac(0, basis.ig_gge0) = std::numeric_limits<double>::quiet_NaN();
    }
    const auto skipped_zero = calculate(false);
    expect_tensor(skipped_zero, actual, 1.0e-8);
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        for (int it = 0; it < cell.ntype; ++it)
        {
            sf.strucFac(it, ig) *= 1.2;
        }
    }
    const auto scaled = calculate(false);
    const double change = scaled(0, 0) - actual(0, 0);
    EXPECT_GT(std::abs(change) * to_kbar(), 1.0);
}

TEST_F(StressEwaTest, GammaWeightOnlyAppliesToPW)
{
    const auto full = calculate(false);
    TestParameters::gamma_only(true);
    const auto lcao_gamma = calculate(false);
    expect_tensor(lcao_gamma, full, 1.0e-8);
    const auto pw_gamma = calculate(true);
    const auto pw_reference = reference_stress(cell, basis, 2.0);
    expect_tensor(pw_gamma, pw_reference, 1.0e-8);

    // Keep one representative of every +/-G pair, including G=0 once.
    TestParameters::gamma_only(false);
    int kept = 0;
    int zero = -1;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const auto n = basis.gdirect[ig];
        const bool positive = n.x > 0 || (n.x == 0 && n.y > 0) || (n.x == 0 && n.y == 0 && n.z >= 0);
        if (positive)
        {
            basis.gdirect[kept] = n;
            basis.gcar[kept] = basis.gcar[ig];
            basis.gg[kept] = basis.gg[ig];
            if (ig == basis.ig_gge0)
            {
                zero = kept;
            }
            ++kept;
        }
    }
    basis.npw = kept;
    basis.ig_gge0 = zero;
    sf.setup(&cell, grid, &basis, false);
    TestParameters::gamma_only(true);
    const auto half_space = calculate(true);
    expect_tensor(half_space, full, 1.0e-8);
}

TEST_F(StressEwaTest, PeriodicTranslationReorderingAndStructureFactorRefresh)
{
    const auto original = calculate(false);
    atoms[0].taud[0].x += 1.0;
    std::reverse(atoms[0].taud.begin(), atoms[0].taud.end());
    update_positions();
    sf.setup(&cell, grid, &basis, false);
    const auto translated = calculate(false);
    expect_tensor(translated, original, 1.0e-8);
    atoms[1].taud[0].y += 0.037;
    const ModuleBase::Matrix3 updated_lattice(1.02, 0.16, 0.04, 0.08, 1.1, 0.17, 0.03, 0.07, 0.95);
    update_geometry(updated_lattice);
    sf.setup(&cell, grid, &basis, false);
    const auto actual = calculate(false);
    const auto expected = reference_stress(cell, basis, 1.0);
    expect_tensor(actual, expected, 1.0e-8);
}

TEST_F(StressEwaTest, SingleSpecies)
{
    cell.ntype = 1;
    cell.nat = 2;
    cell.set_iat2itia();
    sf.setup(&cell, grid, &basis, false);
    const auto actual = calculate(false);
    const auto expected = reference_stress(cell, basis, 1.0);
    expect_tensor(actual, expected, 1.0e-8);
}

TEST_F(StressEwaTest, GZeroAndRealSpaceBelongToTheirOwningProcess)
{
    // G=0 and real space must be counted once even when the owner is not rank 0.
    ModuleBase::matrix previous(3, 3);
    for (const int owner : {0, process_count - 1})
    {
        basis.npw = process_rank == owner ? 1 : 0;
        basis.ig_gge0 = process_rank == owner ? 0 : -1;
        if (basis.npw > 0)
        {
            basis.gdirect[0] = ModuleBase::Vector3<double>(0, 0, 0);
            basis.gcar[0] = ModuleBase::Vector3<double>(0, 0, 0);
            basis.gg[0] = 0;
        }
        sf.setup(&cell, grid, &basis, false);
        const auto actual = calculate(false);
        const auto expected = reference_stress(cell, basis, 1.0);
        expect_tensor(actual, expected, 1.0e-8);
        if (owner == 0)
        {
            previous = actual;
        }
        else
        {
            expect_tensor(actual, previous, 1.0e-8);
        }
    }
}

TEST_F(StressEwaTest, Bspline20AndFractionalGridTranslation)
{
    sf.set(&basis, 20);
    for (int shift = 0; shift < 2; ++shift)
    {
        sf.setup(&cell, grid, &basis, false);
        const auto actual = calculate(false);
        const auto expected = reference_stress(cell, basis, 1.0);
        expect_tensor(actual, expected, 1.0e-4);
        for (auto& atom : atoms)
        {
            for (auto& s : atom.taud)
            {
                s += ModuleBase::Vector3<double>(0.173 / 36, 0.319 / 36, 0.271 / 36);
            }
        }
        update_positions();
    }
}

TEST_F(StressEwaTest, EnergyStrainDerivativeInSixDirections)
{
    const auto lattice = cell.latvec;
    for (const int order : {0, 20})
    {
        sf.set(&basis, order);
        sf.setup(&cell, grid, &basis, false);
        const auto analytic = calculate(false);
        const double volume = cell.omega;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                // Off-diagonal strain has B_ij=B_ji=1/2, so dE/dh=-Omega*sigma_ij.
                for (const double h : {3.0e-4, 1.0e-4})
                {
                    double energies[2] = {};
                    for (int side = 0; side < 2; ++side)
                    {
                        const double step = side == 0 ? -h : h;
                        const ModuleBase::Matrix3 deformed = strained_lattice(lattice, i, j, step);
                        update_geometry(deformed);
                        // Fixed integer G set, grid, fractional positions and cutoff metadata.
                        // SF is invariant under this affine strain and is reused on both sides.
                        std::ofstream unused;
                        energies[side] = H_Ewald_pw::compute_ewald(cell, &basis, sf.strucFac, 0, unused);
                    }
                    const double fd = -(energies[1] - energies[0]) / (2 * h * volume);
                    EXPECT_NEAR(fd * to_kbar(), analytic(i, j) * to_kbar(), 1.0e-3)
                        << "order=" << order << " direction=" << i << j << " h=" << h;
                }
            }
        }
        update_geometry(lattice);
    }
}

TEST_F(StressEwaTest, RebuiltCellBasisAndSplineEnergyStrainDerivative)
{
    const auto lattice = cell.latvec;
    sf.set(&basis, 20);
    sf.setup(&cell, grid, &basis, false);
    const auto analytic = calculate(false);
    const double volume = cell.omega;
    const double h = 1.0e-4;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double energies[2] = {};
            for (int side = 0; side < 2; ++side)
            {
                const double step = side == 0 ? -h : h;
                const ModuleBase::Matrix3 deformed = strained_lattice(lattice, i, j, step);
                update_geometry(deformed);
                // Rebuild the actual FFT distribution and cutoff-selected G set.
                basis.initgrids(cell.lat0, cell.latvec, 80.0);
                basis.initparameters(false, 80.0);
                basis.setuptransform();
                basis.collect_local_pw();
                sf.setup(&cell, grid, &basis, false);
                std::ofstream unused;
                energies[side] = H_Ewald_pw::compute_ewald(cell, &basis, sf.strucFac, 0, unused);
            }
            const double fd = -(energies[1] - energies[0]) / (2 * h * volume);
            EXPECT_NEAR(fd * to_kbar(), analytic(i, j) * to_kbar(), 1.0e-3) << i << "," << j;
        }
    }
}
} // namespace

int main(int argc, char** argv)
{
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, process_count, threads, process_rank);
    int bnd_size = 1;
    int bnd_rank = 0;
    int bnd_group = 0;
    int pool_size = 1;
    int pool_rank = 0;
    int pool = 0;
    Parallel_Global::init_pools(process_count, process_rank, 1, 1,
                               bnd_size, bnd_rank, bnd_group, pool_size, pool_rank, pool);
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
