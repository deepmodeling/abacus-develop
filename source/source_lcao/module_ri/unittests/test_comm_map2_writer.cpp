#include "../ri_2d_comm.h"

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <mpi.h>
#include <set>
#include <utility>

namespace
{
using Cell = std::array<int, 3>;
using Second = std::pair<int, Cell>;
using Block = std::map<int, std::map<Second, double>>;

bool close(const double lhs, const double rhs)
{
    return std::abs(lhs - rhs) < 1.0e-14;
}
} // namespace

// The optional write_cv diagnostic requests the complete atom sets only on
// the writer. This regression owns the contributing blocks on non-writer
// ranks (one of them contributing nothing) and checks that the writer still
// receives every block while the non-writers get empty maps and every input
// map stays unchanged.
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int mpi_rank = 0;
    int mpi_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size != 4)
    {
        if (mpi_rank == 0)
        {
            std::cout << "SKIP: this regression requires four MPI ranks\n";
        }
        MPI_Finalize();
        return 0;
    }

    Block local;
    if (mpi_rank == 1)
    {
        local[10][Second{20, {0, 0, 0}}] = 5.0;
    }
    if (mpi_rank == 2)
    {
        local[11][Second{21, {0, 0, 0}}] = 7.0;
    }
    const Block local_before = local;

    std::set<int> request_atoms0;
    std::set<int> request_atoms1;
    if (mpi_rank == 0)
    {
        request_atoms0 = {10, 11};
        request_atoms1 = {20, 21};
    }
    const Block gathered = RI_2D_Comm::comm_map2_first(MPI_COMM_WORLD, local, request_atoms0, request_atoms1);

    bool passed = local == local_before;
    if (mpi_rank == 0)
    {
        passed = passed && gathered.size() == 2 && gathered.count(10) == 1 && gathered.count(11) == 1;
        if (passed)
        {
            const auto& first_block = gathered.at(10);
            const auto& second_block = gathered.at(11);
            passed = first_block.size() == 1 && second_block.size() == 1
                     && first_block.count(Second{20, {0, 0, 0}}) == 1
                     && second_block.count(Second{21, {0, 0, 0}}) == 1
                     && close(first_block.at(Second{20, {0, 0, 0}}), 5.0)
                     && close(second_block.at(Second{21, {0, 0, 0}}), 7.0);
        }
    }
    else
    {
        passed = passed && gathered.empty();
    }

    int local_ok = passed ? 1 : 0;
    int all_ok = 0;
    MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (mpi_rank == 0)
    {
        std::cout << "writer-only comm_map2_first gather: " << (all_ok == 1 ? "ok" : "FAILED") << '\n';
    }

    MPI_Finalize();
    return all_ok == 1 ? 0 : 1;
}
