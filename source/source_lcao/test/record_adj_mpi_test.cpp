#include "../record_adj.h"

#include "gtest/gtest.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include "source_cell/module_neighbor/test/synthetic_neighbor_unitcell.h"
#include "source_cell/module_neighbor/sltk_atom_arrange.h"
#include "source_cell/read_stru.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <string>
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
using AdjacentRecord = std::array<int, 5>;
using AdjacentRecordList = std::vector<std::vector<AdjacentRecord>>;

int mpi_rank()
{
#ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
#else
    return 0;
#endif
}

int mpi_size()
{
#ifdef __MPI
    int size = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
#else
    return 1;
#endif
}

AdjacentRecordList collect_records(const Record_adj& record_adj)
{
    AdjacentRecordList records(record_adj.na_proc);
    for (int iat = 0; iat < record_adj.na_proc; ++iat)
    {
        for (int adjacent = 0; adjacent < record_adj.na_each[iat]; ++adjacent)
        {
            records[iat].push_back(AdjacentRecord{{record_adj.info[iat][adjacent][0],
                                                   record_adj.info[iat][adjacent][1],
                                                   record_adj.info[iat][adjacent][2],
                                                   record_adj.info[iat][adjacent][3],
                                                   record_adj.info[iat][adjacent][4]}});
        }
    }
    return records;
}

void expect_parallel_metadata(const UnitCell& ucell,
                              const Record_adj& serial_record,
                              const Parallel_Orbitals& parallel_orbitals)
{
    int expected_nnr = 0;
    for (int iat = 0; iat < ucell.nat; ++iat)
    {
        EXPECT_EQ(parallel_orbitals.nlocstart[iat], expected_nnr);
        int expected_dimension = 0;
        const int type1 = ucell.iat2it[iat];
        const int natom1 = ucell.iat2ia[iat];
        const int start1 = ucell.itiaiw2iwt(type1, natom1, 0);
        for (int adjacent = 0; adjacent < serial_record.na_each[iat]; ++adjacent)
        {
            const int type2 = serial_record.info[iat][adjacent][3];
            const int natom2 = serial_record.info[iat][adjacent][4];
            const int start2 = ucell.itiaiw2iwt(type2, natom2, 0);
            for (int orbital1 = 0; orbital1 < ucell.atoms[type1].nw; ++orbital1)
            {
                if (parallel_orbitals.global2local_row(start1 + orbital1) < 0)
                {
                    continue;
                }
                for (int orbital2 = 0; orbital2 < ucell.atoms[type2].nw; ++orbital2)
                {
                    if (parallel_orbitals.global2local_col(start2 + orbital2) >= 0)
                    {
                        ++expected_dimension;
                    }
                }
            }
        }
        EXPECT_EQ(parallel_orbitals.nlocdim[iat], expected_dimension);
        expected_nnr += expected_dimension;
    }
    EXPECT_EQ(parallel_orbitals.nnr, expected_nnr);
}

void expect_mpi_record_matches_serial(SyntheticNeighborCase test_case)
{
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);
    ucell->set_iat2itia();
    int max_atom_count = 0;
    for (int type = 0; type < ucell->ntype; ++type)
    {
        max_atom_count = std::max(max_atom_count, ucell->atoms[type].na);
    }
    ucell->itia2iat.create(ucell->ntype, max_atom_count);
    for (int iat = 0; iat < ucell->nat; ++iat)
    {
        ucell->itia2iat(ucell->iat2it[iat], ucell->iat2ia[iat]) = iat;
    }
    ucell->set_iat2iwt(1);

    const int nlocal = ucell->get_iat2iwt()[ucell->nat - 1]
                       + ucell->atoms[ucell->iat2it[ucell->nat - 1]].nw;
    const double search_radius = test_case.radii.front();
    const double skin_radius = 0.4;
    const double build_radius = search_radius + skin_radius;
    const std::vector<double> orbital_cutoff(ucell->ntype, search_radius * ucell->lat0);
    const std::string output_name
        = "record_adj_mpi_" + test_case.name + "_" + std::to_string(mpi_rank()) + ".out";
    std::ofstream output(output_name);
    Grid_Driver reference_grid(0, 0);
    reference_grid.init(output, *ucell, build_radius, true);
    reference_grid.set_query_radius(search_radius);
    output.close();

    const std::string mpi_output_name
        = "record_adj_prebuilt_mpi_" + test_case.name + "_" + std::to_string(mpi_rank()) + ".out";
    std::ofstream mpi_output(mpi_output_name);
    Grid_Driver prebuilt_mpi_grid(0, 0);
    ModuleNeighbor::MpiGhostExchangeStats stats;
#ifdef __MPI
    atom_arrange::search_mpi(true,
                             mpi_output,
                             prebuilt_mpi_grid,
                             *ucell,
                             search_radius * ucell->lat0,
                             MPI_COMM_WORLD,
                             &stats,
                             skin_radius * ucell->lat0);
#else
    atom_arrange::search_mpi(true,
                             mpi_output,
                             prebuilt_mpi_grid,
                             *ucell,
                             search_radius * ucell->lat0,
                             ModuleNeighbor::kSerialMpiCommWorld,
                             &stats,
                             skin_radius * ucell->lat0);
#endif
    mpi_output.close();

    const std::string direct_output_name
        = "record_adj_direct_mpi_" + test_case.name + "_" + std::to_string(mpi_rank()) + ".out";
    std::ofstream direct_output(direct_output_name);
    Grid_Driver direct_mpi_grid(0, 0);
#ifdef __MPI
    direct_mpi_grid.init_mpi(direct_output, *ucell, build_radius, true, MPI_COMM_WORLD);
#else
    direct_mpi_grid.init_mpi(direct_output,
                             *ucell,
                             build_radius,
                             true,
                             ModuleNeighbor::kSerialMpiCommWorld);
#endif
    direct_mpi_grid.set_query_radius(search_radius);
    direct_output.close();
    EXPECT_EQ(prebuilt_mpi_grid.neighbor_pairs, direct_mpi_grid.neighbor_pairs);
    EXPECT_EQ(prebuilt_mpi_grid.local_center_mask, direct_mpi_grid.local_center_mask);

    Parallel_Orbitals serial_orbitals;
    serial_orbitals.set_serial(nlocal, nlocal);
    serial_orbitals.set_atomic_trace(ucell->get_iat2iwt(), ucell->nat, nlocal);
    Record_adj serial_record;
    serial_record.for_2d(*ucell, reference_grid, serial_orbitals, false, orbital_cutoff);

    Parallel_Orbitals mpi_orbitals;
#ifdef __MPI
    mpi_orbitals.init(nlocal, nlocal, 1, MPI_COMM_WORLD);
#else
    mpi_orbitals.set_serial(nlocal, nlocal);
#endif
    mpi_orbitals.set_atomic_trace(ucell->get_iat2iwt(), ucell->nat, nlocal);
    Record_adj mpi_record;
    mpi_record.for_2d(*ucell, reference_grid, mpi_orbitals, false, orbital_cutoff);

    EXPECT_EQ(mpi_record.na_proc, ucell->nat);
    EXPECT_EQ(collect_records(mpi_record), collect_records(serial_record));
    expect_parallel_metadata(*ucell, serial_record, mpi_orbitals);

    Record_adj prebuilt_record;
    prebuilt_record.for_2d(*ucell,
                           reference_grid,
                           mpi_orbitals,
                           false,
                           orbital_cutoff,
                           &prebuilt_mpi_grid);
    EXPECT_EQ(collect_records(prebuilt_record), collect_records(serial_record));
    expect_parallel_metadata(*ucell, serial_record, mpi_orbitals);

    // A second build verifies that the previous pointer grid is released before
    // the globally reconstructed MPI records replace it.
    prebuilt_record.for_2d(*ucell,
                           reference_grid,
                           mpi_orbitals,
                           false,
                           orbital_cutoff,
                           &prebuilt_mpi_grid);
    EXPECT_EQ(collect_records(prebuilt_record), collect_records(serial_record));
    expect_parallel_metadata(*ucell, serial_record, mpi_orbitals);

    std::remove(output_name.c_str());
    std::remove(mpi_output_name.c_str());
    std::remove(direct_output_name.c_str());
    delete ucell;
}

void expect_prebuilt_grid_validation()
{
    SyntheticNeighborCase test_case = make_synthetic_neighbor_cases()[0];
    UnitCell* ucell = test_case.prepare.SetUcellInfo();
    unitcell::check_dtau(ucell->atoms, ucell->ntype, ucell->lat0, ucell->latvec);
    ucell->set_iat2itia();
    int max_atom_count = 0;
    for (int type = 0; type < ucell->ntype; ++type)
    {
        max_atom_count = std::max(max_atom_count, ucell->atoms[type].na);
    }
    ucell->itia2iat.create(ucell->ntype, max_atom_count);
    for (int iat = 0; iat < ucell->nat; ++iat)
    {
        ucell->itia2iat(ucell->iat2it[iat], ucell->iat2ia[iat]) = iat;
    }
    ucell->set_iat2iwt(1);

    const int nlocal = ucell->get_iat2iwt()[ucell->nat - 1]
                       + ucell->atoms[ucell->iat2it[ucell->nat - 1]].nw;
    const double radius = test_case.radii[0];
    const double build_radius = radius + 0.4;
    const std::vector<double> orbital_cutoff(ucell->ntype, radius * ucell->lat0);
    std::ofstream output("record_adj_validation_" + std::to_string(mpi_rank()) + ".out");
    Grid_Driver reference_grid(0, 0);
    reference_grid.init(output, *ucell, build_radius, true);
    reference_grid.set_query_radius(radius);

    Parallel_Orbitals mpi_orbitals;
#ifdef __MPI
    mpi_orbitals.init(nlocal, nlocal, 1, MPI_COMM_WORLD);
#else
    mpi_orbitals.set_serial(nlocal, nlocal);
#endif
    mpi_orbitals.set_atomic_trace(ucell->get_iat2iwt(), ucell->nat, nlocal);

    Grid_Driver wrong_radius_grid(0, 0);
    Grid_Driver wrong_query_radius_grid(0, 0);
    Grid_Driver wrong_pbc_grid(0, 0);
#ifdef __MPI
    wrong_radius_grid.init_mpi(output, *ucell, build_radius + 0.1, true, MPI_COMM_WORLD);
    wrong_query_radius_grid.init_mpi(output, *ucell, build_radius, true, MPI_COMM_WORLD);
    wrong_pbc_grid.init_mpi(output, *ucell, build_radius, false, MPI_COMM_WORLD);
#else
    wrong_radius_grid.init_mpi(output,
                               *ucell,
                               build_radius + 0.1,
                               true,
                               ModuleNeighbor::kSerialMpiCommWorld);
    wrong_query_radius_grid.init_mpi(output,
                                     *ucell,
                                     build_radius,
                                     true,
                                     ModuleNeighbor::kSerialMpiCommWorld);
    wrong_pbc_grid.init_mpi(output,
                            *ucell,
                            build_radius,
                            false,
                            ModuleNeighbor::kSerialMpiCommWorld);
#endif
    wrong_radius_grid.set_query_radius(radius);
    wrong_query_radius_grid.set_query_radius(radius + 0.1);
    wrong_pbc_grid.set_query_radius(radius);
    Record_adj record;
    EXPECT_THROW(record.for_2d(*ucell,
                               reference_grid,
                               mpi_orbitals,
                               false,
                               orbital_cutoff,
                               &wrong_radius_grid),
                 std::invalid_argument);
    EXPECT_THROW(record.for_2d(*ucell,
                               reference_grid,
                               mpi_orbitals,
                               false,
                               orbital_cutoff,
                               &wrong_query_radius_grid),
                 std::invalid_argument);
    EXPECT_THROW(record.for_2d(*ucell,
                               reference_grid,
                               mpi_orbitals,
                               false,
                               orbital_cutoff,
                               &wrong_pbc_grid),
                 std::invalid_argument);

#ifdef __MPI
    if (mpi_size() > 1)
    {
        Grid_Driver wrong_communicator_grid(0, 0);
        wrong_communicator_grid.init_mpi(output, *ucell, build_radius, true, MPI_COMM_SELF);
        wrong_communicator_grid.set_query_radius(radius);
        EXPECT_THROW(record.for_2d(*ucell,
                                   reference_grid,
                                   mpi_orbitals,
                                   false,
                                   orbital_cutoff,
                                   &wrong_communicator_grid),
                     std::invalid_argument);
    }
#endif

    output.close();
    std::remove(("record_adj_validation_" + std::to_string(mpi_rank()) + ".out").c_str());
    delete ucell;
}
} // namespace

TEST(RecordAdjMpiTest, OrthogonalPeriodicMatchesSerial)
{
    expect_mpi_record_matches_serial(make_synthetic_neighbor_cases()[0]);
}

TEST(RecordAdjMpiTest, TriclinicPeriodicMatchesSerial)
{
    expect_mpi_record_matches_serial(make_synthetic_neighbor_cases()[1]);
}

TEST(RecordAdjMpiTest, RejectsMismatchedPrebuiltGrid)
{
    expect_prebuilt_grid_validation();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
#ifdef __MPI
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(&argc, &argv);
    }
#endif
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized)
    {
        MPI_Finalize();
    }
#endif
    return result;
}
