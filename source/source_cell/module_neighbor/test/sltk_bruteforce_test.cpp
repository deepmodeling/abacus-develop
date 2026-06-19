#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "../sltk_grid_driver.h"
#include "prepare_unitcell.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <tuple>
#include <valarray>
#include <vector>

#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
LCAO_Orbitals::LCAO_Orbitals()
{
}
LCAO_Orbitals::~LCAO_Orbitals()
{
}
#endif

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
    this->start_mag = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_mag;
}

namespace
{
struct NeighborKey
{
    int type = 0;
    int atom = 0;
    int box_x = 0;
    int box_y = 0;
    int box_z = 0;

    bool operator<(const NeighborKey& other) const
    {
        return std::tie(type, atom, box_x, box_y, box_z)
               < std::tie(other.type, other.atom, other.box_x, other.box_y, other.box_z);
    }

    bool operator==(const NeighborKey& other) const
    {
        return std::tie(type, atom, box_x, box_y, box_z)
               == std::tie(other.type, other.atom, other.box_x, other.box_y, other.box_z);
    }
};

struct ExpandLayers
{
    int plus_x = 0;
    int minus_x = 0;
    int plus_y = 0;
    int minus_y = 0;
    int plus_z = 0;
    int minus_z = 0;
};

ModuleBase::Matrix3 MatrixFromRows(const std::vector<double>& rows)
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = rows[0];
    latvec.e12 = rows[1];
    latvec.e13 = rows[2];
    latvec.e21 = rows[3];
    latvec.e22 = rows[4];
    latvec.e23 = rows[5];
    latvec.e31 = rows[6];
    latvec.e32 = rows[7];
    latvec.e33 = rows[8];
    return latvec;
}

std::valarray<double> ValarrayFromVector(const std::vector<double>& values)
{
    return std::valarray<double>(values.data(), values.size());
}

std::valarray<int> ValarrayFromVector(const std::vector<int>& values)
{
    return std::valarray<int>(values.data(), values.size());
}

std::unique_ptr<UnitCell> MakeCell(const double lat0,
                                   const ModuleBase::Matrix3& latvec,
                                   const std::vector<std::string>& labels,
                                   const std::vector<int>& atoms_per_type,
                                   const std::vector<double>& coords)
{
    std::vector<std::string> pp_files;
    std::vector<std::string> pp_types;
    std::vector<std::string> orb_files;
    std::vector<double> masses;
    for (const std::string& label: labels)
    {
        pp_files.push_back(label + ".upf");
        pp_types.push_back("upf201");
        orb_files.push_back(label + ".orb");
        masses.push_back(1.0);
    }

    const std::vector<double> latvec_rows = {latvec.e11,
                                             latvec.e12,
                                             latvec.e13,
                                             latvec.e21,
                                             latvec.e22,
                                             latvec.e23,
                                             latvec.e31,
                                             latvec.e32,
                                             latvec.e33};
    UcellTestPrepare prep("custom",
                          2,
                          false,
                          false,
                          false,
                          "volume",
                          lat0,
                          ValarrayFromVector(latvec_rows),
                          labels,
                          pp_files,
                          pp_types,
                          orb_files,
                          ValarrayFromVector(atoms_per_type),
                          masses,
                          "Cartesian",
                          ValarrayFromVector(coords));
    return std::unique_ptr<UnitCell>(prep.SetUcellInfo());
}

double CrossNorm(const ModuleBase::Vector3<double>& lhs, const ModuleBase::Vector3<double>& rhs)
{
    const double x = lhs.y * rhs.z - lhs.z * rhs.y;
    const double y = lhs.z * rhs.x - lhs.x * rhs.z;
    const double z = lhs.x * rhs.y - lhs.y * rhs.x;
    return std::sqrt(x * x + y * y + z * z);
}

ExpandLayers BruteForceLayers(const UnitCell& ucell, const double radius_lat0)
{
    const ModuleBase::Vector3<double> a1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    const ModuleBase::Vector3<double> a2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    const ModuleBase::Vector3<double> a3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);
    const double lat0_volume = ucell.lat0 * ucell.lat0 * ucell.lat0;

    const int dx = std::ceil(CrossNorm(a2, a3) * radius_lat0 / ucell.omega * lat0_volume);
    const int dy = std::ceil(CrossNorm(a3, a1) * radius_lat0 / ucell.omega * lat0_volume);
    const int dz = std::ceil(CrossNorm(a1, a2) * radius_lat0 / ucell.omega * lat0_volume);

    return {dx + 1, dx, dy + 1, dy, dz + 1, dz};
}

std::vector<NeighborKey> BruteForceNeighbors(const UnitCell& ucell,
                                             const double radius_lat0,
                                             const int center_type,
                                             const int center_atom)
{
    const ExpandLayers layers = BruteForceLayers(ucell, radius_lat0);
    const ModuleBase::Vector3<double> a1(ucell.latvec.e11, ucell.latvec.e12, ucell.latvec.e13);
    const ModuleBase::Vector3<double> a2(ucell.latvec.e21, ucell.latvec.e22, ucell.latvec.e23);
    const ModuleBase::Vector3<double> a3(ucell.latvec.e31, ucell.latvec.e32, ucell.latvec.e33);
    const ModuleBase::Vector3<double>& center = ucell.atoms[center_type].tau[center_atom];
    const double radius2 = radius_lat0 * radius_lat0;

    std::vector<NeighborKey> result;
    for (int ix = -layers.minus_x; ix < layers.plus_x; ++ix)
    {
        for (int iy = -layers.minus_y; iy < layers.plus_y; ++iy)
        {
            for (int iz = -layers.minus_z; iz < layers.plus_z; ++iz)
            {
                const ModuleBase::Vector3<double> shift(ix * a1.x + iy * a2.x + iz * a3.x,
                                                        ix * a1.y + iy * a2.y + iz * a3.y,
                                                        ix * a1.z + iy * a2.z + iz * a3.z);
                for (int it = 0; it < ucell.ntype; ++it)
                {
                    for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
                    {
                        const ModuleBase::Vector3<double> image(ucell.atoms[it].tau[ia].x + shift.x,
                                                                ucell.atoms[it].tau[ia].y + shift.y,
                                                                ucell.atoms[it].tau[ia].z + shift.z);
                        const double dx = center.x - image.x;
                        const double dy = center.y - image.y;
                        const double dz = center.z - image.z;
                        const double dr2 = dx * dx + dy * dy + dz * dz;
                        const bool central_self
                            = (it == center_type && ia == center_atom && ix == 0 && iy == 0 && iz == 0);
                        if (!central_self && dr2 <= radius2)
                        {
                            result.push_back({it, ia, ix, iy, iz});
                        }
                    }
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<NeighborKey> GridNeighbors(const Grid_Driver& grid,
                                       const UnitCell& ucell,
                                       const int center_type,
                                       const int center_atom)
{
    AdjacentAtomInfo adjs;
    grid.Find_atom(ucell, center_type, center_atom, &adjs);
    EXPECT_EQ(adjs.ntype.size(), static_cast<size_t>(adjs.adj_num + 1));
    EXPECT_EQ(adjs.natom.size(), static_cast<size_t>(adjs.adj_num + 1));
    EXPECT_EQ(adjs.box.size(), static_cast<size_t>(adjs.adj_num + 1));
    EXPECT_EQ(adjs.adjacent_tau.size(), static_cast<size_t>(adjs.adj_num + 1));
    EXPECT_EQ(adjs.ntype[adjs.adj_num], center_type);
    EXPECT_EQ(adjs.natom[adjs.adj_num], center_atom);
    EXPECT_EQ(adjs.box[adjs.adj_num].x, 0);
    EXPECT_EQ(adjs.box[adjs.adj_num].y, 0);
    EXPECT_EQ(adjs.box[adjs.adj_num].z, 0);

    std::vector<NeighborKey> result;
    for (int i = 0; i < adjs.adj_num; ++i)
    {
        result.push_back({adjs.ntype[i], adjs.natom[i], adjs.box[i].x, adjs.box[i].y, adjs.box[i].z});
    }
    std::sort(result.begin(), result.end());
    return result;
}

void ExpectGridMatchesBruteForce(const UnitCell& ucell,
                                 const double radius_lat0,
                                 const std::string& output_name)
{
    Grid_Driver grid(false, false);
    std::ofstream ofs(output_name);
    grid.init(ofs, ucell, radius_lat0, true);
    ofs.close();

    const ExpandLayers layers = BruteForceLayers(ucell, radius_lat0);
    EXPECT_EQ(grid.getGlayerX(), layers.plus_x);
    EXPECT_EQ(grid.getGlayerX_minus(), layers.minus_x);
    EXPECT_EQ(grid.getGlayerY(), layers.plus_y);
    EXPECT_EQ(grid.getGlayerY_minus(), layers.minus_y);
    EXPECT_EQ(grid.getGlayerZ(), layers.plus_z);
    EXPECT_EQ(grid.getGlayerZ_minus(), layers.minus_z);

    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            EXPECT_EQ(GridNeighbors(grid, ucell, it, ia), BruteForceNeighbors(ucell, radius_lat0, it, ia))
                << "center type=" << it << " atom=" << ia;
        }
    }
    std::remove(output_name.c_str());
}

bool ContainsNeighbor(const std::vector<NeighborKey>& neighbors, const NeighborKey& key)
{
    return std::find(neighbors.begin(), neighbors.end(), key) != neighbors.end();
}

void ExpectReciprocalGridNeighbors(const UnitCell& ucell,
                                   const double radius_lat0,
                                   const std::string& output_name)
{
    Grid_Driver grid(false, false);
    std::ofstream ofs(output_name);
    grid.init(ofs, ucell, radius_lat0, true);
    ofs.close();

    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            const std::vector<NeighborKey> neighbors = GridNeighbors(grid, ucell, it, ia);
            for (const NeighborKey& neighbor: neighbors)
            {
                const std::vector<NeighborKey> reverse_neighbors = GridNeighbors(grid, ucell, neighbor.type, neighbor.atom);
                const NeighborKey reverse_key{it, ia, -neighbor.box_x, -neighbor.box_y, -neighbor.box_z};
                EXPECT_TRUE(ContainsNeighbor(reverse_neighbors, reverse_key))
                    << "missing reverse neighbor for center type=" << it << " atom=" << ia;
            }
        }
    }
    std::remove(output_name.c_str());
}
} // namespace

TEST(SltkBruteForceTest, OrthogonalMultiTypeCellMatchesIndependentBruteForce)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B"},
                                {2, 2},
                                {0.05, 0.05, 0.05,
                                 0.55, 0.50, 0.50,
                                 0.95, 0.05, 0.05,
                                 0.12, 0.06, 0.05});

    ExpectGridMatchesBruteForce(*ucell, 0.18, "bruteforce_orthogonal.out");
}

TEST(SltkBruteForceTest, PbcFaceEdgeAndCornerImagesMatchIndependentBruteForce)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A"},
                                {4},
                                {0.02, 0.02, 0.02,
                                 0.98, 0.02, 0.02,
                                 0.98, 0.98, 0.02,
                                 0.98, 0.98, 0.98});

    ExpectGridMatchesBruteForce(*ucell, 0.08, "bruteforce_pbc_edge.out");
}

TEST(SltkBruteForceTest, TriclinicCellMatchesIndependentBruteForce)
{
    const auto ucell = MakeCell(25.0,
                                MatrixFromRows({1.0, 0.22, 0.07, 0.15, 1.0, 0.18, 0.08, 0.25, 1.0}),
                                {"A", "B"},
                                {2, 1},
                                {0.03, 0.96, 0.04,
                                 0.52, 0.48, 0.51,
                                 0.94, 0.08, 0.06});

    ExpectGridMatchesBruteForce(*ucell, 0.24, "bruteforce_triclinic.out");
}

TEST(SltkBruteForceTest, CutoffBoundaryIncludesEqualDistanceAndRejectsOutside)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B"},
                                {2, 2},
                                {0.50, 0.50, 0.50,
                                 0.75, 0.50, 0.50,
                                 0.25, 0.50, 0.50,
                                 0.751, 0.50, 0.50});

    ExpectGridMatchesBruteForce(*ucell, 0.25, "bruteforce_cutoff_boundary.out");
}

TEST(SltkBruteForceTest, DenseClusterAcrossPeriodicBoundaryKeepsImageIdentities)
{
    const auto ucell = MakeCell(30.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B", "C"},
                                {3, 2, 2},
                                {0.015, 0.020, 0.025,
                                 0.985, 0.020, 0.025,
                                 0.020, 0.985, 0.025,
                                 0.985, 0.985, 0.025,
                                 0.030, 0.025, 0.985,
                                 0.975, 0.975, 0.975,
                                 0.500, 0.500, 0.500});

    ExpectGridMatchesBruteForce(*ucell, 0.09, "bruteforce_dense_boundary_cluster.out");
}

TEST(SltkBruteForceTest, StronglySkewedTriclinicCellMatchesIndependentBruteForce)
{
    const auto ucell = MakeCell(18.0,
                                MatrixFromRows({1.0, 0.46, 0.19, 0.38, 1.0, 0.31, 0.12, 0.27, 1.0}),
                                {"A", "B"},
                                {3, 3},
                                {0.06, 0.08, 0.91,
                                 0.18, 0.88, 0.12,
                                 0.43, 0.47, 0.51,
                                 0.92, 0.09, 0.08,
                                 0.78, 0.82, 0.87,
                                 0.51, 0.49, 0.48});

    ExpectGridMatchesBruteForce(*ucell, 0.22, "bruteforce_strongly_skewed.out");
}

TEST(SltkBruteForceTest, DeterministicIrregularSmallCellMatchesIndependentBruteForce)
{
    const auto ucell = MakeCell(22.0,
                                MatrixFromRows({1.0, 0.17, 0.03, 0.09, 1.0, 0.21, 0.04, 0.13, 1.0}),
                                {"A", "B", "C"},
                                {2, 3, 2},
                                {0.137, 0.941, 0.219,
                                 0.773, 0.113, 0.617,
                                 0.951, 0.887, 0.071,
                                 0.064, 0.481, 0.902,
                                 0.382, 0.284, 0.336,
                                 0.621, 0.742, 0.158,
                                 0.229, 0.059, 0.809});

    ExpectGridMatchesBruteForce(*ucell, 0.19, "bruteforce_irregular_small.out");
}

TEST(SltkBruteForceTest, SameAtomPeriodicImagesAcrossMultipleShellsAreKept)
{
    const auto ucell = MakeCell(12.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A"},
                                {1},
                                {0.50, 0.50, 0.50});

    ExpectGridMatchesBruteForce(*ucell, 1.05, "bruteforce_same_atom_images.out");
}

TEST(SltkBruteForceTest, AnisotropicCellUsesDifferentExpansionLayersPerDirection)
{
    const auto ucell = MakeCell(16.0,
                                MatrixFromRows({0.55, 0.00, 0.00, 0.00, 1.70, 0.00, 0.00, 0.00, 1.25}),
                                {"A", "B"},
                                {2, 2},
                                {0.08, 0.10, 0.12,
                                 0.91, 0.11, 0.13,
                                 0.49, 1.48, 0.60,
                                 0.52, 0.24, 1.10});

    ExpectGridMatchesBruteForce(*ucell, 0.34, "bruteforce_anisotropic.out");
}

TEST(SltkBruteForceTest, NearCutoffAcrossPeriodicFaceKeepsOnlyInsideImage)
{
    const auto ucell = MakeCell(10.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A"},
                                {4},
                                {0.0100, 0.5000, 0.5000,
                                 0.9899, 0.5000, 0.5000,
                                 0.6500, 0.5000, 0.5000,
                                 0.6502, 0.5000, 0.5000});

    ExpectGridMatchesBruteForce(*ucell, 0.35, "bruteforce_periodic_cutoff_face.out");
}

TEST(SltkBruteForceTest, ManyTypesWithUnevenAtomCountsPreserveTypeAndAtomIdentity)
{
    const auto ucell = MakeCell(24.0,
                                MatrixFromRows({1.0, 0.11, 0.02, 0.05, 1.0, 0.16, 0.03, 0.08, 1.0}),
                                {"A", "B", "C", "D"},
                                {1, 3, 1, 2},
                                {0.04, 0.05, 0.06,
                                 0.96, 0.05, 0.06,
                                 0.50, 0.50, 0.50,
                                 0.21, 0.78, 0.32,
                                 0.08, 0.94, 0.07,
                                 0.92, 0.91, 0.93,
                                 0.43, 0.18, 0.88});

    ExpectGridMatchesBruteForce(*ucell, 0.17, "bruteforce_many_types_uneven.out");
}

TEST(SltkBruteForceTest, VacuumLikeSparseCellDoesNotInventNeighbors)
{
    const auto ucell = MakeCell(80.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B"},
                                {2, 1},
                                {0.10, 0.10, 0.10,
                                 0.80, 0.80, 0.80,
                                 0.45, 0.45, 0.45});

    ExpectGridMatchesBruteForce(*ucell, 0.05, "bruteforce_vacuum_sparse.out");
}

TEST(SltkBruteForceTest, LeadingEmptyAtomTypeDoesNotBreakIndexing)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.08, 0.02, 0.03, 1.0, 0.06, 0.01, 0.04, 1.0}),
                                {"Empty", "A", "B"},
                                {0, 2, 1},
                                {0.05, 0.05, 0.05,
                                 0.95, 0.05, 0.05,
                                 0.52, 0.50, 0.50});

    ExpectGridMatchesBruteForce(*ucell, 0.16, "bruteforce_leading_empty_type.out");
}

TEST(SltkBruteForceTest, VerySmallNonzeroSeparationIsNotMistakenForSelf)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B"},
                                {1, 1},
                                {0.40, 0.40, 0.40,
                                 0.400001, 0.40, 0.40});

    ExpectGridMatchesBruteForce(*ucell, 0.00001, "bruteforce_tiny_separation.out");
}

TEST(SltkBruteForceTest, SubAngstromScaleSeparationsRemainDistinctDownToDoublePrecision)
{
    const std::vector<double> separations = {1.0e-8, 1.0e-10, 1.0e-12, 1.0e-14, 1.0e-15};
    for (const double delta: separations)
    {
        const auto ucell = MakeCell(20.0,
                                    MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                    {"A", "B"},
                                    {1, 1},
                                    {0.40, 0.40, 0.40,
                                     0.40 + delta, 0.40, 0.40});

        ExpectGridMatchesBruteForce(*ucell, delta * 2.0, "bruteforce_double_precision_separation.out");
    }
}

TEST(SltkBruteForceTest, DirectionalTinySeparationsInSkewedCellRemainDistinct)
{
    const std::vector<double> separations = {1.0e-8, 1.0e-10, 1.0e-12, 1.0e-14};
    for (const double delta: separations)
    {
        const auto ucell = MakeCell(20.0,
                                    MatrixFromRows({1.0, 0.31, 0.07, 0.18, 1.0, 0.23, 0.05, 0.29, 1.0}),
                                    {"A"},
                                    {8},
                                    {0.40, 0.40, 0.40,
                                     0.40 + delta, 0.40, 0.40,
                                     0.40, 0.40 + delta, 0.40,
                                     0.40, 0.40, 0.40 + delta,
                                     0.40 + delta, 0.40 + delta, 0.40,
                                     0.40 + delta, 0.40, 0.40 + delta,
                                     0.40, 0.40 + delta, 0.40 + delta,
                                     0.40 + delta, 0.40 + delta, 0.40 + delta});

        ExpectGridMatchesBruteForce(*ucell, delta * 3.0, "bruteforce_directional_tiny_separation.out");
    }
}

TEST(SltkBruteForceTest, DeterministicFuzzSmallCellsMatchIndependentBruteForce)
{
    const std::vector<ModuleBase::Matrix3> lattices = {
        MatrixFromRows({1.0, 0.00, 0.00, 0.00, 1.0, 0.00, 0.00, 0.00, 1.0}),
        MatrixFromRows({1.0, 0.18, 0.04, 0.07, 1.0, 0.11, 0.03, 0.09, 1.0}),
        MatrixFromRows({0.82, 0.28, 0.05, 0.16, 1.13, 0.18, 0.04, 0.22, 0.91}),
    };
    const std::vector<double> radii = {0.13, 0.19, 0.27};

    for (size_t ilat = 0; ilat < lattices.size(); ++ilat)
    {
        for (size_t ir = 0; ir < radii.size(); ++ir)
        {
            std::vector<double> coords;
            coords.reserve(18);
            for (int i = 0; i < 6; ++i)
            {
                const double seed = static_cast<double>((ilat + 1) * 37 + (ir + 1) * 19 + i * 23);
                coords.push_back(std::fmod(0.137 * seed + 0.071 * i, 0.92) + 0.04);
                coords.push_back(std::fmod(0.193 * seed + 0.053 * i, 0.92) + 0.04);
                coords.push_back(std::fmod(0.271 * seed + 0.029 * i, 0.92) + 0.04);
            }

            const auto ucell = MakeCell(24.0,
                                        lattices[ilat],
                                        {"A", "B", "C"},
                                        {2, 3, 1},
                                        coords);
            ExpectGridMatchesBruteForce(*ucell, radii[ir], "bruteforce_deterministic_fuzz.out");
        }
    }
}

TEST(SltkBruteForceTest, NearZeroAndOneFractionalBoundariesMatchIndependentBruteForce)
{
    const double eps = 1.0e-12;
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.16, 0.03, 0.05, 1.0, 0.12, 0.02, 0.07, 1.0}),
                                {"A", "B"},
                                {4, 4},
                                {eps, 0.50, 0.50,
                                 1.0 - eps, 0.50, 0.50,
                                 0.50, eps, 0.50,
                                 0.50, 1.0 - eps, 0.50,
                                 0.50, 0.50, eps,
                                 0.50, 0.50, 1.0 - eps,
                                 eps, eps, eps,
                                 1.0 - eps, 1.0 - eps, 1.0 - eps});

    ExpectGridMatchesBruteForce(*ucell, 0.09, "bruteforce_near_boundary_normalization.out");
}

TEST(SltkBruteForceTest, ThinStronglySkewedCellMatchesIndependentBruteForce)
{
    const auto ucell = MakeCell(18.0,
                                MatrixFromRows({1.0, 0.65, 0.00, 0.15, 1.0, 0.00, 0.03, 0.04, 0.22}),
                                {"A", "B"},
                                {3, 3},
                                {0.05, 0.06, 0.08,
                                 0.95, 0.07, 0.09,
                                 0.48, 0.52, 0.11,
                                 0.08, 0.91, 0.17,
                                 0.89, 0.86, 0.19,
                                 0.51, 0.49, 0.03});

    ExpectGridMatchesBruteForce(*ucell, 0.11, "bruteforce_thin_strongly_skewed.out");
}

TEST(SltkBruteForceTest, NeighborRelationIsReciprocalWithOppositePeriodicBox)
{
    const auto ucell = MakeCell(18.0,
                                MatrixFromRows({1.0, 0.25, 0.05, 0.12, 1.0, 0.17, 0.04, 0.20, 1.0}),
                                {"A", "B", "C"},
                                {2, 2, 1},
                                {0.04, 0.05, 0.06,
                                 0.96, 0.07, 0.06,
                                 0.10, 0.91, 0.08,
                                 0.88, 0.86, 0.92,
                                 0.52, 0.48, 0.50});

    ExpectReciprocalGridNeighbors(*ucell, 0.19, "bruteforce_reciprocal.out");
}

TEST(SltkBruteForceTest, DISABLED_PeriodicBoundaryEquivalentPositionsShouldRemainDistinctNeighbors)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B"},
                                {1, 1},
                                {0.0, 0.50, 0.50,
                                 1.0, 0.50, 0.50});

    ExpectGridMatchesBruteForce(*ucell, 0.10, "bruteforce_periodic_equivalent_positions.out");
}

TEST(SltkBruteForceTest, DISABLED_DistinctAtomsAtSamePositionShouldBeReportedAsNeighbors)
{
    const auto ucell = MakeCell(20.0,
                                MatrixFromRows({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}),
                                {"A", "B"},
                                {1, 1},
                                {0.40, 0.40, 0.40,
                                 0.40, 0.40, 0.40});

    ExpectGridMatchesBruteForce(*ucell, 0.10, "bruteforce_overlapping_atoms.out");
}
