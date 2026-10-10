#include "../symm_rotation.h"
#include "gtest/gtest.h"

// K-point generation is outside this test: use explicit stars, but provide
// the virtual symbols needed by the existing lightweight rotation test target.
void ModuleCell::ReciprocalGrid::renew(const int&)
{
    ADD_FAILURE() << "Unexpected k-point generation";
}
void K_Vectors::renew(const int&)
{
    ADD_FAILURE() << "Unexpected k-point generation";
}
void K_Vectors::reduce_by_symmetry(const UnitCell&, const ModuleSymmetry::Symmetry&,
                                  bool, std::string&, bool&, const int, std::ofstream&)
{
    ADD_FAILURE() << "Unexpected k-point reduction";
}

namespace
{
using Complex = std::complex<double>;

// Independent dense product in the stored (transposed density) convention.
std::vector<Complex> rotate_reference(const std::vector<Complex>& density,
                                      const std::vector<Complex>& rotation,
                                      const int n)
{
    std::vector<Complex> result(n * n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int a = 0; a < n; ++a)
            {
                for (int b = 0; b < n; ++b)
                {
                    result[i + j * n] += rotation[a + i * n] * density[a + b * n]
                                         * std::conj(rotation[b + j * n]);
                }
            }
        }
    }
    return result;
}

void check_little_group_restoration(const int nspin)
{
    const int channels = nspin == 2 ? 2 : 1;
    const int n = 4;
    Parallel_2D pv;
    pv.init(n, n, 1, MPI_COMM_WORLD);
    ModuleSymmetry::Symmetry_rotation rotation;
    std::vector<Complex> identity(n * n, 0.0);
    std::vector<Complex> little(n * n, 0.0);
    std::vector<Complex> representative(n * n, 0.0);
    std::vector<Complex> alternate(n * n, 0.0);
    const int sign[n] = {1, 1, -1, -1};
    for (int i = 0; i < n; ++i)
    {
        identity[i + i * n] = 1.0;
        little[i + i * n] = sign[i];
        const int row = (i + 1) % n;
        representative[row + i * n] = std::polar(1.0, 0.3 * i);
        alternate[row + i * n] = double(sign[row]) * representative[row + i * n];
    }
    auto local = [&pv, n](const std::vector<Complex>& dense) {
        std::vector<Complex> result(pv.get_local_size());
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if (pv.in_this_processor(i, j))
                {
                    result[pv.global2local_row(i) + pv.global2local_col(j) * pv.get_row_size()]
                        = dense[i + j * n];
                }
            }
        }
        return result;
    };
    rotation.set_density_rotations_for_testing(
        {{{0, local(identity)}, {1, local(little)}, {2, local(representative)}, {3, local(alternate)}}},
        {{0, 1}}, 4, nspin);
    K_Vectors kv;
    kv.set_nkstot(channels);
    kv.set_nkstot_nospin(2);
    // Single-pool (KPAR=1) scenario: kv.get_nks() (local) equals kv.get_nkstot() (global),
    // and ik2iktot is the identity map. With only one global ibz-k here, any value mod
    // kv.kstars.size()==1 is 0, so the exact ik2iktot values don't matter, only its size.
    kv.set_nks(channels);
    kv.ik2iktot.assign(channels, 0);
    kv.kstars = {{{0, {0.25, 0.0, 0.0}}, {2, {0.0, 0.25, 0.0}}}};
    std::vector<std::vector<Complex>> inputs;
    std::vector<std::vector<Complex>> expected;
    for (int spin = 0; spin < channels; ++spin)
    {
        std::vector<Complex> density(n * n);
        std::vector<Complex> projected(n * n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                const Complex value((spin + 1) * (2.0 + i + j), 0.2 * (i - j));
                density[i + j * n] = value;
                // For this C2 little group, averaging removes exactly the odd blocks.
                projected[i + j * n] = sign[i] == sign[j] ? 0.5 * value : Complex(0.0);
            }
        }
        inputs.push_back(local(density));
        expected.push_back(local(projected));
        expected.push_back(local(rotate_reference(projected, representative, n)));
    }
    const auto restored = rotation.restore_dm(kv, inputs, pv);
    ASSERT_EQ(restored.size(), expected.size());
    for (size_t k = 0; k < expected.size(); ++k)
    {
        for (size_t i = 0; i < expected[k].size(); ++i)
        {
            EXPECT_NEAR(std::abs(restored[k][i] - expected[k][i]), 0.0, 1e-12);
        }
    }
    // A different representative of the same star must give the same density.
    kv.kstars = {{{0, {0.25, 0.0, 0.0}}, {3, {0.0, 0.25, 0.0}}}};
    const auto changed_representative = rotation.restore_dm(kv, inputs, pv);
    for (size_t k = 0; k < expected.size(); ++k)
    {
        for (size_t i = 0; i < expected[k].size(); ++i)
        {
            EXPECT_NEAR(std::abs(changed_representative[k][i] - restored[k][i]), 0.0, 1e-12);
        }
    }
}
} // namespace

TEST(SymmetryDensityRestoration, LittleGroupAndStarWeight)
{
    check_little_group_restoration(1);
}

TEST(SymmetryDensityRestoration, IndependentSpinChannels)
{
    check_little_group_restoration(2);
}

TEST(SymmetryDensityRestoration, SpinorDensity)
{
    check_little_group_restoration(4);
}

namespace
{
class CellSymmetryRotation : public ModuleSymmetry::Symmetry_rotation
{
  public:
    const std::vector<std::map<int, std::vector<Complex>>>& rotations() const
    {
        return this->Ms_;
    }
};
}

TEST(SymmetryDensityRestoration, RebuildAfterCellSymmetryAnalysis)
{
    // A common translation crosses the cell boundary for only part of the basis.
    // The point group stays unchanged, but the lattice returns and Bloch phases change.
    UnitCell cell;
    std::vector<Atom> atoms(2);
    cell.atoms = atoms.data();
    cell.ntype = 2;
    cell.nat = 3;
    cell.lmax = 0;
    cell.latvec = ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    cell.a1 = ModuleBase::Vector3<double>(1, 0, 0);
    cell.a2 = ModuleBase::Vector3<double>(0, 1, 0);
    cell.a3 = ModuleBase::Vector3<double>(0, 0, 1);
    cell.st.iat2it = new int[3]{0, 0, 1};
    cell.st.iat2ia = new int[3]{0, 1, 0};
    atoms[0].label = "A";
    atoms[0].na = 2;
    atoms[0].nw = 1;
    atoms[0].iw2l = {0};
    atoms[0].stapos_wf = 0;
    atoms[0].taud = {{0.1, 0.2, 0.3}, {0.4, 0.3, 0.2}};
    atoms[0].tau = atoms[0].taud;
    atoms[1].label = "B";
    atoms[1].na = 1;
    atoms[1].nw = 1;
    atoms[1].iw2l = {0};
    atoms[1].stapos_wf = 2;
    atoms[1].taud = {{0.25, 0.25, 0.25}};
    atoms[1].tau = atoms[1].taud;
    std::ofstream log;
    const int representation[2] = {0, 0};
    const std::string calculation = "cell-relax";
    cell.symm.analy_sys(cell.lat, cell.st, cell.atoms, log, 1e-6, 1, calculation, representation);
    const int old_operations = cell.symm.nrotk;
    K_Vectors kv;
    kv.set_nks(1);
    kv.set_nkstot(1);
    kv.set_nkstot_nospin(2);
    int total_kpoints = 1;
    kv.para_k.kinfo(total_kpoints, 1, 0, 0, 1, 1);
    kv.ik2iktot = {0};
    kv.kvec_d = {{0.25, 0, 0}};
    kv.kstars = {{{0, {0.25, 0, 0}}}};
    for (int operation = 0; operation < cell.symm.nrotk; ++operation)
    {
        kv.kstars[0][operation] = kv.kvec_d[0] * cell.symm.kgmatrix[operation];
    }
    Parallel_2D pv;
    pv.init(3, 3, 1, MPI_COMM_WORLD);
    const ModuleSymmetry::TC period = {2, 2, 2};
    const std::vector<ModuleSymmetry::TC> cells = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CellSymmetryRotation reused;
    reused.find_irreducible_sector(cell.symm, cell.atoms, cell.st, cells, period, cell.lat);
    reused.cal_Ms(kv, cell, pv, 1);
    const auto old_rotations = reused.rotations();

    atoms[0].taud = {{0.9, 0.2, 0.3}, {0.2, 0.3, 0.2}};
    atoms[1].taud = {{0.05, 0.25, 0.25}};
    atoms[0].tau = atoms[0].taud;
    atoms[1].tau = atoms[1].taud;
    cell.symm.analy_sys(cell.lat, cell.st, cell.atoms, log, 1e-6, 1, calculation, representation);
    ASSERT_EQ(cell.symm.nrotk, old_operations);
    reused.reset_symmetry();
    reused.find_irreducible_sector(cell.symm, cell.atoms, cell.st, cells, period, cell.lat);
    reused.cal_Ms(kv, cell, pv, 1);
    CellSymmetryRotation fresh;
    fresh.find_irreducible_sector(cell.symm, cell.atoms, cell.st, cells, period, cell.lat);
    fresh.cal_Ms(kv, cell, pv, 1);
    EXPECT_EQ(reused.get_irreducible_sector(), fresh.get_irreducible_sector());
    EXPECT_EQ(reused.rotations(), fresh.rotations());
    EXPECT_NE(old_rotations, fresh.rotations());
    for (int atom = 0; atom < cell.nat; ++atom)
    {
        for (int operation = 0; operation < cell.symm.nrotk; ++operation)
        {
            const auto actual = reused.get_return_lattice(atom, operation);
            const auto expected = fresh.get_return_lattice(atom, operation);
            EXPECT_EQ(actual.x, expected.x);
            EXPECT_EQ(actual.y, expected.y);
            EXPECT_EQ(actual.z, expected.z);
        }
    }
}

TEST(SymmetryDensityRestoration, ReturnLatticePreservesBoundaryRepresentatives)
{
    ModuleSymmetry::Symmetry symmetry;
    symmetry.epsilon = 3.2e-5;
    ModuleSymmetry::Irreducible_Sector sector;
    const ModuleBase::Matrix3 reflection(-1, 0, 0, 0, 1, 0, 0, 0, 1);
    const ModuleBase::Vector3<double> translation(0.0, 0.0, 0.0);
    const ModuleBase::Vector3<double> source(0.99999, 0.25, 0.5);
    const ModuleBase::Vector3<double> mapped(0.00001, 0.25, 0.5);
    const auto lattice = sector.get_return_lattice(symmetry, reflection, translation, source, mapped);
    EXPECT_DOUBLE_EQ(lattice.x, -1.0);
    EXPECT_DOUBLE_EQ(lattice.y, 0.0);
    EXPECT_DOUBLE_EQ(lattice.z, 0.0);

    // Changing either atom's representative must change its integer lattice shift.
    const ModuleBase::Vector3<double> shifted_source(1.99999, 0.25, 0.5);
    const ModuleBase::Vector3<double> shifted_mapped(1.00001, 0.25, 0.5);
    const auto source_lattice = sector.get_return_lattice(
        symmetry, reflection, translation, shifted_source, mapped);
    const auto mapped_lattice = sector.get_return_lattice(
        symmetry, reflection, translation, source, shifted_mapped);
    EXPECT_DOUBLE_EQ(source_lattice.x, -2.0);
    EXPECT_DOUBLE_EQ(mapped_lattice.x, -2.0);
}

TEST(SymmetryDensityRestoration, ReturnLatticePreservesTranslationRepresentative)
{
    ModuleSymmetry::Symmetry symmetry;
    symmetry.epsilon = 3.2e-5;
    ModuleSymmetry::Irreducible_Sector sector;
    const ModuleBase::Matrix3 identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
    const ModuleBase::Vector3<double> source(0.25, 0.5, 0.75);
    const ModuleBase::Vector3<double> translation(0.99999, 0.0, 0.0);
    const ModuleBase::Vector3<double> mapped(0.24999, 0.5, 0.75);
    const auto lattice = sector.get_return_lattice(symmetry, identity, translation, source, mapped);
    EXPECT_DOUBLE_EQ(lattice.x, 1.0);
    EXPECT_DOUBLE_EQ(lattice.y, 0.0);
    EXPECT_DOUBLE_EQ(lattice.z, 0.0);
}
