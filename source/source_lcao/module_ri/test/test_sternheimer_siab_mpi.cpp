#include "source_lcao/module_ri/sternheimer_chi0_mpi.h"
#include "source_lcao/module_ri/sternheimer_siab_mpi.h"
#include "source_lcao/module_ri/sternheimer_siab_overlap.h"
#include "source_lcao/module_ri/sternheimer_siab_writer.h"
#include "source_lcao/module_ri/sternheimer_rpa.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <iterator>
#include <stdexcept>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace
{

namespace siab = module_ri::sternheimer_siab;
using Complex = std::complex<double>;

std::string read_text(const std::string& path)
{
    std::ifstream input(path.c_str(), std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

siab::Provenance test_provenance()
{
    siab::Provenance provenance;
    provenance.abacus_commit = "19ab21e01d02cc805604ed77a6e269af698fdd1d";
    provenance.auxiliary_basis_sha256 = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    provenance.cell_bohr = {8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0};
    provenance.ecut_ry = 25.0;
    provenance.kernel = "full_coulomb";
    provenance.orbital_sha256 = "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";
    provenance.pseudopotential_sha256 = "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc";
    provenance.spin_convention = "occupation_in_metadata";
    return provenance;
}

siab::ReferenceRow make_row(const int frequency_index,
                            const double frequency_ha,
                            const std::vector<Complex>& y,
                            const std::vector<std::vector<Complex>>& primitives)
{
    siab::ReferenceRow row;
    row.occupied_state = 0;
    row.auxiliary_channel = 0;
    row.frequency_index = frequency_index;
    row.frequency_ha = frequency_ha;
    row.occupation = 2.0;
    row.frequency_weight = 0.25 + 0.5 * frequency_index;
    row.norm = siab::norm(y, 0.5);
    row.q = siab::overlap_q(y, primitives, 0.5);
    return row;
}

siab::ReferenceRow make_global_row(const int occupied_state,
                                   const int frequency_index,
                                   const int auxiliary_channel)
{
    const double tag = 100.0 * occupied_state + 10.0 * frequency_index + auxiliary_channel;
    siab::ReferenceRow row;
    row.occupied_state = occupied_state;
    row.auxiliary_channel = auxiliary_channel;
    row.frequency_index = frequency_index;
    row.frequency_ha = 0.2 + 0.6 * frequency_index;
    row.occupation = 2.0;
    row.frequency_weight = 0.25 + 0.5 * frequency_index;
    row.norm = 1.0 + tag;
    row.q = {Complex(tag + 1.0, -tag), Complex(0.5 * tag, tag + 2.0)};
    return row;
}

TEST(SternheimerSIABMPI, ReassemblesPWSlabsInFullGridOrder)
{
    siab::PrimitiveSlab rank0;
    rank0.startz = 0;
    rank0.nplane = 1;
    rank0.values = {{Complex(1.0, 0.5), Complex(3.0, 1.5)}, {Complex(10.0, -0.5), Complex(30.0, -1.5)}};

    siab::PrimitiveSlab rank1;
    rank1.startz = 1;
    rank1.nplane = 1;
    rank1.values = {{Complex(2.0, 1.0), Complex(4.0, 2.0)}, {Complex(20.0, -1.0), Complex(40.0, -2.0)}};

    const auto full = siab::assemble_full_primitive_grids({rank1, rank0}, 2, 2);
    ASSERT_EQ(full.size(), 2);
    EXPECT_EQ(full[0], (std::vector<Complex>{{1.0, 0.5}, {2.0, 1.0}, {3.0, 1.5}, {4.0, 2.0}}));
    EXPECT_EQ(full[1], (std::vector<Complex>{{10.0, -0.5}, {20.0, -1.0}, {30.0, -1.5}, {40.0, -2.0}}));
}

TEST(SternheimerSIABMPI, RejectsDuplicateOrMissingPWPlanes)
{
    siab::PrimitiveSlab first;
    first.startz = 0;
    first.nplane = 1;
    first.values = {{Complex(1.0, 0.0)}};

    siab::PrimitiveSlab duplicate = first;
    EXPECT_THROW(siab::assemble_full_primitive_grids({first, duplicate}, 1, 2), std::invalid_argument);
    EXPECT_THROW(siab::assemble_full_primitive_grids({first}, 1, 2), std::invalid_argument);
}

#ifdef __MPI
TEST(SternheimerSIABMPI, FullPrimitiveAndVariableRowsMatchSerialWithoutDoubleCounting)
{
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2)
    {
        GTEST_SKIP() << "This regression requires exactly two MPI ranks.";
    }

    // PW local ordering is [ixy][local iz].  Each rank owns one z plane.
    const int startz = rank;
    const int nplane = 1;
    std::vector<std::vector<Complex>> local_primitives(2, std::vector<Complex>(2));
    local_primitives[0][0] = Complex(1.0 + rank, 0.5 * (rank + 1));
    local_primitives[0][1] = Complex(3.0 + rank, 1.5 + 0.5 * rank);
    local_primitives[1][0] = Complex(10.0 * (rank + 1), -0.5 * (rank + 1));
    local_primitives[1][1] = Complex(30.0 + 10.0 * rank, -1.5 - 0.5 * rank);

    const auto full_primitives
        = siab::allgather_full_primitive_grids(local_primitives, 2, 2, startz, nplane, MPI_COMM_WORLD);
    ASSERT_EQ(full_primitives.size(), 2);
    ASSERT_EQ(full_primitives[0].size(), 4);

    const std::vector<Complex> y0{{1.0, 1.0}, {2.0, -1.0}, {0.5, 0.25}, {-1.0, 0.75}};
    const std::vector<Complex> y1{{-0.5, 0.5}, {1.5, 0.25}, {2.0, -0.5}, {0.25, 1.0}};

    // Deliberately reverse owner/frequency order.  Y is replicated and complete;
    // only the owning rank forms its row, with no overlap Allreduce.
    std::vector<siab::ReferenceRow> local_rows;
    if (rank == 0)
    {
        local_rows.push_back(make_row(1, 0.8, y1, full_primitives));
    }
    else
    {
        local_rows.push_back(make_row(0, 0.2, y0, full_primitives));
    }

    const auto gathered = siab::gather_reference_rows_to_root(local_rows, 2, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        ASSERT_EQ(gathered.size(), 2);
        EXPECT_EQ(gathered[0].frequency_index, 0);
        EXPECT_EQ(gathered[1].frequency_index, 1);

        const auto expected0 = make_row(0, 0.2, y0, full_primitives);
        const auto expected1 = make_row(1, 0.8, y1, full_primitives);
        EXPECT_NEAR(gathered[0].norm, expected0.norm, 1.0e-12);
        EXPECT_NEAR(gathered[1].norm, expected1.norm, 1.0e-12);
        for (std::size_t ie = 0; ie != 2; ++ie)
        {
            EXPECT_NEAR(std::abs(gathered[0].q[ie] - expected0.q[ie]), 0.0, 1.0e-12);
            EXPECT_NEAR(std::abs(gathered[1].q[ie] - expected1.q[ie]), 0.0, 1.0e-12);
        }

        const std::vector<siab::PrimitiveBlock> blocks{{"H", 0, 0, 0, 2, 0}};
        const std::vector<Complex> overlap_s = siab::overlap_s(full_primitives, 0.5);
        const std::string mpi_path = "sternheimer_siab_mpi_rows.dat";
        const std::string serial_path = "sternheimer_siab_serial_rows.dat";
        siab::write_v1(mpi_path, 0.5, blocks, gathered, overlap_s, test_provenance());
        siab::write_v1(serial_path, 0.5, blocks, {expected1, expected0}, overlap_s, test_provenance());
        EXPECT_EQ(read_text(mpi_path), read_text(serial_path));
        std::remove(mpi_path.c_str());
        std::remove(serial_path.c_str());
    }
    else
    {
        EXPECT_TRUE(gathered.empty());
    }
}

TEST(SternheimerSIABMPI, SupportsRankWithNoLocalReferenceRows)
{
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2)
    {
        GTEST_SKIP() << "This regression requires exactly two MPI ranks.";
    }

    std::vector<siab::ReferenceRow> local_rows;
    if (rank == 1)
    {
        const std::vector<std::vector<Complex>> primitive{{Complex(1.0, 0.0)}};
        local_rows.push_back(make_row(1, 0.8, {Complex(2.0, 0.0)}, primitive));
        local_rows.push_back(make_row(0, 0.2, {Complex(1.0, 0.0)}, primitive));
    }

    const auto gathered = siab::gather_reference_rows_to_root(local_rows, 1, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        ASSERT_EQ(gathered.size(), 2);
        EXPECT_EQ(gathered[0].frequency_index, 0);
        EXPECT_EQ(gathered[1].frequency_index, 1);
    }
}

TEST(SternheimerSIABMPI, GlobalEquationOwnershipWritesCanonicalSerialFile)
{
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2)
    {
        GTEST_SKIP() << "This regression requires exactly two MPI ranks.";
    }

    constexpr int occupied_count = 2;
    constexpr int frequency_count = 2;
    constexpr int channel_count = 4;
    std::vector<siab::ReferenceRow> local_rows;
    std::vector<siab::ReferenceRow> serial_rows;
    for (int occupied = occupied_count - 1; occupied >= 0; --occupied)
    {
        for (int frequency = frequency_count - 1; frequency >= 0; --frequency)
        {
            for (int channel = channel_count - 1; channel >= 0; --channel)
            {
                const siab::ReferenceRow row = make_global_row(occupied, frequency, channel);
                serial_rows.push_back(row);
                const int owner = ModuleRI::SternheimerRPA::global_equation_owner(occupied,
                                                                                  frequency,
                                                                                  channel,
                                                                                  frequency_count,
                                                                                  channel_count,
                                                                                  size,
                                                                                  1);
                if (owner == rank)
                {
                    local_rows.push_back(row);
                }
            }
        }
    }
    ASSERT_EQ(local_rows.size(), 8);

    const auto gathered = siab::gather_reference_rows_to_root(local_rows, 2, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        ASSERT_EQ(gathered.size(), 16);
        const std::vector<siab::PrimitiveBlock> blocks{{"H", 0, 0, 0, 2, 0}};
        const std::vector<Complex> overlap_s{Complex(1.0, 0.0),
                                             Complex(0.0, 0.0),
                                             Complex(0.0, 0.0),
                                             Complex(1.0, 0.0)};
        const std::string mpi_path = "sternheimer_siab_global_mpi_rows.dat";
        const std::string serial_path = "sternheimer_siab_global_serial_rows.dat";
        siab::write_v1(mpi_path, 0.5, blocks, gathered, overlap_s, test_provenance());
        siab::write_v1(serial_path, 0.5, blocks, serial_rows, overlap_s, test_provenance());
        EXPECT_EQ(read_text(mpi_path), read_text(serial_path));
        std::remove(mpi_path.c_str());
        std::remove(serial_path.c_str());
    }
    else
    {
        EXPECT_TRUE(gathered.empty());
    }
}

TEST(SternheimerSIABMPI, GlobalEquationChi0ReductionMatchesSerialBranch)
{
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2)
    {
        GTEST_SKIP() << "This regression requires exactly two MPI ranks.";
    }

    std::vector<Complex> local_branch(4, Complex(0.0, 0.0));
    if (rank == 0)
    {
        local_branch[0] = Complex(1.0, 0.5);
        local_branch[2] = Complex(3.0, -2.0);
    }
    else
    {
        local_branch[1] = Complex(2.0, 1.0);
        local_branch[3] = Complex(4.0, -0.25);
    }

    module_ri::sternheimer_chi0::reduce_branch_to_root(local_branch, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        const std::vector<Complex> serial_branch{Complex(1.0, 0.5),
                                                  Complex(2.0, 1.0),
                                                  Complex(3.0, -2.0),
                                                  Complex(4.0, -0.25)};
        EXPECT_EQ(local_branch, serial_branch);
        EXPECT_EQ(ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(local_branch, 2),
                  ModuleRI::SternheimerRPA::symmetrize_chi0_imaginary_frequency(serial_branch, 2));
    }
    else
    {
        EXPECT_TRUE(local_branch.empty());
    }
}

#endif

} // namespace

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
