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

namespace
{
template <typename T>
std::vector<T> transform_for_test(const std::vector<T>& old, const std::vector<double>& ao,
    const std::vector<T>& current, int naos, int bands);
template <typename T>
std::vector<T> project(const std::vector<T>& previous, const std::vector<T>& overlap, int no, int nv);
}

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
                if (shifted.norm() < 1.6) { expected.push_back(shifted.norm()); }
            }
        }
    }
    std::sort(actual.begin(), actual.end());
    std::sort(expected.begin(), expected.end());
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) { EXPECT_NEAR(actual[i], expected[i], 1e-13); }
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
    const auto mo = transform_for_test(coefficients, ao, coefficients, 2, 2);
    EXPECT_EQ(mo, ao);
}

TEST(RootOverlap, MOTransformConjugatesOldOrbitals)
{
    using C = std::complex<double>;
    const std::vector<C> old{C(0, 1), C(0), C(0), C(1)};
    const std::vector<C> current{C(1), C(0), C(0), C(0, 1)};
    const std::vector<double> ao{1, 0.2, -0.4, 0.9};
    const auto mo = transform_for_test(old, ao, current, 2, 2);
    EXPECT_EQ(mo[0], C(0, -1));
    EXPECT_EQ(mo[1], C(0.2));
    EXPECT_EQ(mo[2], C(-0.4));
    EXPECT_EQ(mo[3], C(0, 0.9));
}


namespace
{
void setup(Parallel_2D& target, int rows, int cols, const Parallel_2D& grid)
{
#ifdef __MPI
    target.set(rows, cols, 1, grid.blacs_ctxt);
#else
    target.set_serial(rows, cols);
#endif
}

template <typename T>
std::vector<T> scatter(const std::vector<T>& full, const Parallel_2D& layout)
{
    std::vector<T> local(layout.get_local_size(), T(0));
    const int columns = layout.get_global_col_size();
    for (int col = 0; col < layout.get_col_size(); ++col)
    {
        const int global_col = layout.local2global_col(col);
        for (int row = 0; row < layout.get_row_size(); ++row)
        {
            const int global_row = layout.local2global_row(row);
            local[col * layout.get_row_size() + row] = full[global_row * columns + global_col];
        }
    }
    return local;
}

// Independent test oracles: direct formula contractions, with no BLAS/PBLAS calls.
template <typename T> T conjugate(T value) { return value; }
template <> std::complex<double> conjugate(std::complex<double> value) { return std::conj(value); }

template <typename T>
std::vector<T> reference_mo(const std::vector<T>& old, const std::vector<double>& ao,
    const std::vector<T>& current, int naos, int bands)
{
    std::vector<T> result(bands * bands, T(0));
    for (int i = 0; i < bands; ++i)
    {
        for (int j = 0; j < bands; ++j)
        {
            for (int mu = 0; mu < naos; ++mu)
            {
                for (int nu = 0; nu < naos; ++nu)
                {
                    result[i * bands + j] += conjugate(old[i * naos + mu])
                        * ao[mu * naos + nu] * current[j * naos + nu];
                }
            }
        }
    }
    return result;
}

template <typename T>
std::vector<T> reference_projection(const std::vector<T>& previous,
    const std::vector<T>& overlap, int no, int nv)
{
    const int bands = no + nv;
    std::vector<T> result(no * nv, T(0));
    for (int j = 0; j < no; ++j)
    {
        for (int b = 0; b < nv; ++b)
        {
            for (int i = 0; i < no; ++i)
            {
                for (int a = 0; a < nv; ++a)
                {
                    result[j * nv + b] += overlap[i * bands + j] * previous[i * nv + a]
                        * conjugate(overlap[(no + a) * bands + no + b]);
                }
            }
        }
    }
    return result;
}

template <typename T>
std::vector<T> gather(const std::vector<T>& local, const Parallel_2D& layout)
{
    const int columns = layout.get_global_col_size();
    const int count = layout.get_global_row_size() * columns;
    std::vector<T> result(count, T(0));
    EXPECT_EQ(local.size(), layout.get_local_size());
    for (int col = 0; col < layout.get_col_size(); ++col)
    {
        const int global_col = layout.local2global_col(col);
        for (int row = 0; row < layout.get_row_size(); ++row)
        {
            const int global_row = layout.local2global_row(row);
            result[global_row * columns + global_col] = local[col * layout.get_row_size() + row];
        }
    }
    Parallel_Reduce::reduce_all(result.data(), count);
    return result;
}

template <typename T>
std::vector<T> transform_for_test(const std::vector<T>& old, const std::vector<double>& ao,
    const std::vector<T>& current, int naos, int bands)
{
    Parallel_2D pc;
#ifdef __MPI
    pc.init(naos, bands, 1, MPI_COMM_WORLD);
#else
    pc.set_serial(naos, bands);
#endif
    Parallel_2D ps;
    Parallel_2D po;
    setup(ps, naos, naos, pc);
    setup(po, bands, bands, pc);
    std::vector<T> old_rows(old.size());
    std::vector<T> current_rows(current.size());
    for (int i = 0; i < naos; ++i)
    {
        for (int band = 0; band < bands; ++band)
        {
            old_rows[i * bands + band] = old[band * naos + i];
            current_rows[i * bands + band] = current[band * naos + i];
        }
    }
    const std::vector<T> typed_ao(ao.begin(), ao.end());
    const auto local_old = scatter(old_rows, pc);
    const auto local_current = scatter(current_rows, pc);
    const auto local_ao = scatter(typed_ao, ps);
    const auto actual = LR::mo_overlap_dist(local_old, local_ao, local_current, pc, ps, po);
    return gather(actual, po);
}

template <typename T>
std::vector<T> project(const std::vector<T>& previous, const std::vector<T>& overlap, int no, int nv)
{
    const int bands = no + nv;
    Parallel_2D po;
#ifdef __MPI
    po.init(bands, bands, 1, MPI_COMM_WORLD);
#else
    po.set_serial(bands, bands);
#endif
    Parallel_2D px;
    setup(px, nv, no, po);
    std::vector<T> previous_rows(previous.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a) { previous_rows[a * no + i] = previous[i * nv + a]; }
    }
    const auto local_previous = scatter(previous_rows, px);
    const auto local_overlap = scatter(overlap, po);
    const auto actual = LR::project_reference_dist(local_previous, local_overlap, po, px, no, nv);
    const auto rows = gather(actual, px);
    std::vector<T> result(previous.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a) { result[i * nv + a] = rows[a * no + i]; }
    }
    return result;
}

template <typename T>
void check_mo(int naos, int bands, T phase)
{
    std::vector<double> ao(naos * naos);
    std::vector<T> old(naos * bands);
    std::vector<T> current(naos * bands);
    for (int i = 0; i < naos; ++i)
    {
        for (int j = 0; j < naos; ++j) { ao[i * naos + j] = (2 * i - j + 1) / 7.0; }
        for (int band = 0; band < bands; ++band)
        {
            old[band * naos + i] = phase * ((i + band + 1) / 11.0);
            current[band * naos + i] = T((i - 2 * band + 2) / 13.0);
        }
    }
    const auto actual = transform_for_test(old, ao, current, naos, bands);
    const auto expected = reference_mo(old, ao, current, naos, bands);
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - expected[i]), 0.0, 1e-13);
    }
}
template <typename T>
void check_projection(int no, int nv, T phase)
{
    const int bands = no + nv;
    Parallel_2D po;
#ifdef __MPI
    po.init(bands, bands, 1, MPI_COMM_WORLD);
#else
    po.set_serial(bands, bands);
#endif
    Parallel_2D px;
    setup(px, nv, no, po);
    std::vector<T> overlap(bands * bands);
    std::vector<T> previous(no * nv);
    for (int i = 0; i < bands; ++i)
    {
        for (int j = 0; j < bands; ++j)
        {
            overlap[i * bands + j] = T((i - 2 * j + 1) / 17.0) + phase * ((i + j + 1) / 19.0);
        }
    }
    for (int i = 0; i < no * nv; ++i) { previous[i] = phase * ((i + 1) / 23.0); }
    // Convert occupied-major X into row-major A=X^T before scattering into px.
    std::vector<T> previous_rows(previous.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a) { previous_rows[a * no + i] = previous[i * nv + a]; }
    }
    const auto local_previous = scatter(previous_rows, px);
    const auto local_overlap = scatter(overlap, po);
    const auto actual = LR::project_reference_dist(local_previous, local_overlap, po, px, no, nv);
    const auto expected = reference_projection(previous, overlap, no, nv);
    std::vector<T> expected_rows(expected.size());
    for (int i = 0; i < no; ++i)
    {
        for (int a = 0; a < nv; ++a) { expected_rows[a * no + i] = expected[i * nv + a]; }
    }
    const auto local_expected = scatter(expected_rows, px);
    ASSERT_EQ(actual.size(), local_expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        EXPECT_NEAR(std::abs(actual[i] - local_expected[i]), 0.0, 1e-13);
    }
}

}

TEST(RootDistributed, RealNonsymmetricRectangularTransformation) { check_mo(5, 3, 1.0); }
TEST(RootDistributed, ComplexAdjointTransformation)
{
    const std::complex<double> phase(0.6, 0.8);
    check_mo(5, 3, phase);
}
TEST(RootDistributed, EmptyCoefficientBlocksParticipate) { check_mo(2, 1, 1.0); }

TEST(RootDistributed, RealProjectionPreservesWindowLoss) { check_projection(2, 3, 0.4); }
TEST(RootDistributed, ComplexProjectionUsesVirtualAdjoint)
{
    const std::complex<double> phase(0.3, 0.4);
    check_projection(2, 3, phase);
}
TEST(RootDistributed, EmptyReferenceBlocksParticipate) { check_projection(1, 1, 0.4); }

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
