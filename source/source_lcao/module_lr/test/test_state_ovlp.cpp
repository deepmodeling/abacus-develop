#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_basis/module_nao/two_center_integrator.h"
#include "source_cell/unitcell.h"
#include <memory>
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include <gtest/gtest.h>
#include <complex>
#include "source_lcao/module_lr/state_ovlp.h"
#include "source_base/matrix3.h"
#include <algorithm>
#include <cmath>

#include "state_ovlp_test.h"

TEST(RootOverlap, IncludesPairsAbsentFromBothSameGeometryLists)
{
    // Old centers 0, 1.1; new centers 0.2, 1.3. Both intra-step distances
    // exceed 1, but old second/new first are only 0.9 apart.
    const ModuleBase::Matrix3 lattice(10, 0, 0, 0, 10, 0, 0, 0, 10);
    const ModuleBase::Vector3<double> displacement(-0.9, 0, 0);
    const auto images = LR::cross_image_indices(displacement, lattice, 1.0);
    ASSERT_EQ(images.size(), 1u);
    const auto shifted = displacement + images[0] * lattice;
    EXPECT_NEAR(shifted.x, -0.9, 1e-14);
}

TEST(RootOverlap, IncludesImagesAfterLargeTranslationsInSkewCell)
{
    const ModuleBase::Matrix3 lattice(2, 0, 0, 1.7, 1.1, 0, 0.4, 0.3, 2.2);
    const ModuleBase::Vector3<int> translation(4, -3, 2);
    const ModuleBase::Vector3<double> residual(0.2, 0.1, 0.3);
    const ModuleBase::Vector3<double> displacement = translation * lattice + residual;
    const auto images = LR::cross_image_indices(displacement, lattice, 1.6);
    std::vector<double> actual;
    std::vector<double> expected;
    for (const auto& image : images)
    {
        const auto shifted = displacement + image * lattice;
        const double distance = shifted.norm();
        actual.push_back(distance);
    }
    for (int x = -8; x <= 8; ++x)
    {
        for (int y = -8; y <= 8; ++y)
        {
            for (int z = -8; z <= 8; ++z)
            {
                const ModuleBase::Vector3<int> image(x, y, z);
                const auto shifted = displacement + image * lattice;
                if (shifted.norm() < 1.6)
                {
                    expected.push_back(shifted.norm());
                }
            }
        }
    }
    std::sort(actual.begin(), actual.end());
    std::sort(expected.begin(), expected.end());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_NEAR(actual[i], expected[i], 1e-13);
    }
}

TEST(RootOverlap, RetainsIntegerTranslationsForFourierPhases)
{
    const ModuleBase::Matrix3 lattice(2, 0, 0, 0, 3, 0, 0, 0, 4);
    const ModuleBase::Vector3<double> displacement(2.2, 0.0, 0.0);
    const auto indices = LR::cross_image_indices(displacement, lattice, 0.5);
    ASSERT_EQ(indices.size(), 1u);
    EXPECT_EQ(indices[0].x, -1);
    EXPECT_EQ(indices[0].y, 0);
    EXPECT_EQ(indices[0].z, 0);
    const auto shifted = displacement + indices[0] * lattice;
    EXPECT_NEAR(shifted.x, 0.2, 1e-14);
}

TEST(RootOverlap, RemovesOrbitalSignGaugeFromRootSelection)
{
    const double a = 1.0 / std::sqrt(2.0);
    const std::vector<double> previous{a, a};
    const std::vector<double> overlap{1, 0, 0, 0, 1, 0, 0, 0, -1};
    const auto projected = project(previous, overlap, 1, 2);
    EXPECT_DOUBLE_EQ(projected[0], a);
    EXPECT_DOUBLE_EQ(projected[1], -a);
    const double physical_root = projected[0] * a - projected[1] * a;
    const double other_root = projected[0] * a + projected[1] * a;
    EXPECT_NEAR(physical_root, 1.0, 1e-14);
    EXPECT_NEAR(other_root, 0.0, 1e-14);
}

TEST(RootOverlap, PreservesComplexConjugationAndBothOrbitalRotations)
{
    using C = std::complex<double>;
    const std::vector<C> previous{C(1), C(0), C(0), C(0)};
    // Swap occupied orbitals with a phase and independently rotate virtuals.
    const double a = 1.0 / std::sqrt(2.0);
    const std::vector<C> overlap{C(0), C(0, 1), C(0), C(0),
                                C(1), C(0), C(0), C(0),
                                C(0), C(0), C(a), C(0, a),
                                C(0), C(0), C(0, a), C(a)};
    const auto projected = project(previous, overlap, 2, 2);
    EXPECT_EQ(projected[0], C(0));
    EXPECT_EQ(projected[1], C(0));
    EXPECT_EQ(projected[2], C(0, a));
    EXPECT_EQ(projected[3], C(a));
}

TEST(RootOverlap, RetainsProjectionLossInsteadOfRenormalizing)
{
    const std::vector<double> previous{1};
    const std::vector<double> overlap{0.5, 0.7, -0.3, 0.8};
    const auto projected = project(previous, overlap, 1, 1);
    EXPECT_DOUBLE_EQ(projected[0], 0.4);
}


TEST(RootOverlap, UsesOrderedNonsymmetricAOOverlap)
{
    const std::vector<double> coefficients{1, 0, 0, 1};
    const std::vector<double> ao{1, 0.2, -0.4, 0.9};
    const auto mo = transform_for_test(coefficients, ao, coefficients, 2, 2, 0, 2, 1);
    EXPECT_EQ(mo, ao);
}

TEST(RootOverlap, MOTransformConjugatesOldOrbitals)
{
    using C = std::complex<double>;
    const std::vector<C> old{C(0, 1), C(0), C(0), C(1)};
    const std::vector<C> current{C(1), C(0), C(0), C(0, 1)};
    const std::vector<double> ao{1, 0.2, -0.4, 0.9};
    const auto mo = transform_for_test(old, ao, current, 2, 2, 0, 2, 1);
    EXPECT_EQ(mo[0], C(0, -1));
    EXPECT_EQ(mo[1], C(0.2));
    EXPECT_EQ(mo[2], C(-0.4));
    EXPECT_EQ(mo[3], C(0, 0.9));
}


TEST(RootDistributed, RealNonsymmetricRectangularTransformation)
{
    check_mo(5, 3, 1.0, 0, 3, 1);
}
TEST(RootDistributed, ComplexAdjointTransformation)
{
    const std::complex<double> phase(0.6, 0.8);
    check_mo(5, 3, phase, 0, 3, 1);
}
TEST(RootDistributed, EmptyCoefficientBlocksParticipate)
{
    check_mo(2, 1, 1.0, 0, 1, 1);
}

TEST(RootDistributed, AlignedOccupiedBlock)
{
    check_mo(9, 8, 1.0, 0, 3, 2);
}
TEST(RootDistributed, MisalignedVirtualBlock)
{
    const std::complex<double> phase(0.6, 0.8);
    check_mo(9, 8, phase, 3, 5, 32);
}

TEST(RootDistributed, RealProjectionPreservesWindowLoss)
{
    check_projection(2, 3, 0.4, 1);
}
TEST(RootDistributed, ComplexProjectionUsesVirtualAdjoint)
{
    const std::complex<double> phase(0.3, 0.4);
    check_projection(2, 3, phase, 1);
}
TEST(RootDistributed, EmptyReferenceBlocksParticipate)
{
    check_projection(1, 1, 0.4, 1);
}

TEST(RootDistributed, RealProjectionWithMisalignedBlockTwo)
{
    check_projection(3, 5, 0.4, 2);
}
TEST(RootDistributed, ComplexProjectionWithMisalignedBlockFour)
{
    const std::complex<double> phase(0.3, 0.4);
    check_projection(5, 6, phase, 4);
}
TEST(RootDistributed, ProjectionWithBlockThirtyTwoAndMultipleBlocks)
{
    check_projection(35, 37, 0.4, 32);
}
TEST(RootDistributed, ComplexLargeBlocksWithEmptyRanks)
{
    const std::complex<double> phase(0.3, 0.4);
    check_projection(3, 5, phase, 32);
}

TEST(RootOverlap, CrossGeometryAOOverlapIsOrderedAndNonsymmetric)
{
    UnitCell cell;
    std::unique_ptr<Atom> atoms(new Atom);
    cell.atoms = atoms.get();
    cell.ntype = 1;
    cell.nat = 2;
    cell.lat0 = 1.0;
    cell.latvec = ModuleBase::Matrix3(10, 0, 0, 0, 10, 0, 0, 0, 10);
    cell.iat2it = {0, 0};
    cell.iat2ia = {0, 1};
    atoms->nw = 1;
    atoms->na = 2;
    atoms->iw2l = {0};
    atoms->iw2n = {0};
    atoms->iw2m = {0};
    atoms->tau = {{0.2, 0.0, 0.0}, {1.3, 0.0, 0.0}};
    const std::vector<int> offsets{0, 1};
    cell.set_iat2iwt_for_test(offsets, 1);
    Parallel_Orbitals distribution;
    distribution.set_serial(2, 2);
    distribution.set_atomic_trace(cell.get_iat2iwt(), cell.nat, 2);
    TwoCenterIntegrator integrator;
    const std::vector<ModuleBase::Vector3<double>> old_positions{{0.0, 0.0, 0.0}, {1.1, 0.0, 0.0}};
    const std::vector<double> cutoff{0.5};
    const auto actual = LR::cross_ao_overlap(cell, old_positions, cutoff, integrator, distribution);
    // Existing integral mock returns one: old second/new first is inside the cutoff,
    // whereas old first/new second is outside. Output is column-major.
    const std::vector<double> expected{1.0, 1.0, 0.0, 1.0};
    EXPECT_EQ(actual, expected);
}

int main(int argc, char** argv)
{
    int processes = 1;
    int threads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, processes, threads, rank);
#ifdef __MPI
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
