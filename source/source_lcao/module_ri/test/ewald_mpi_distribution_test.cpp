#include "../ewald_mpi_utils.h"

#include <RI/distribute/Distribute_Equally.h>
#include <RI/distribute/Divide_Atoms.h>
#include <array>
#include <iostream>
#include <mpi.h>
#include <set>
#include <tuple>
#include <vector>

namespace
{
using Cell = std::array<int, 3>;
using Key = std::tuple<int, int, int, int, int>;

template <typename Distribution>
std::set<Key> make_local_keys(const Distribution& distribution)
{
    std::set<Key> keys;
    for (const int atom0: distribution.first)
    {
        for (const auto& atom1_cell: distribution.second.at(0))
        {
            const Cell& cell = atom1_cell.second;
            keys.emplace(atom0, atom1_cell.first, cell[0], cell[1], cell[2]);
        }
    }
    return keys;
}

std::set<Key> make_full_keys(const std::vector<int>& atoms, const Cell& period)
{
    std::set<Key> keys;
    const auto atom_cells = RI::Divide_Atoms::traversal_atom_period(atoms, period);
    for (const int atom0: atoms)
    {
        for (const auto& atom1_cell: atom_cells)
        {
            const Cell& cell = atom1_cell.second;
            keys.emplace(atom0, atom1_cell.first, cell[0], cell[1], cell[2]);
        }
    }
    return keys;
}

std::size_t count_intersection(const std::set<Key>& lhs, const std::set<Key>& rhs)
{
    std::size_t count = 0;
    for (const Key& key: lhs)
    {
        count += rhs.count(key);
    }
    return count;
}
} // namespace

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

    const std::vector<int> atoms = {0, 1, 2};
    const Cell bare_period = {13, 13, 3};
    const Cell gaussian_period = {7, 7, 3};
    const std::set<Key> gaussian_support = make_full_keys(atoms, gaussian_period);

    const auto common_distribution
        = EwaldVqDetail::distribute_common_realspace_tasks(MPI_COMM_WORLD, atoms, bare_period);
    const std::set<Key> local_bare = make_local_keys(common_distribution);

    std::set<Key> local_gaussian_on_bare_owner;
    for (const Key& key: local_bare)
    {
        if (gaussian_support.count(key) != 0)
        {
            local_gaussian_on_bare_owner.insert(key);
        }
    }

    const auto independently_distributed_gaussian
        = RI::Distribute_Equally::distribute_atoms_periods(MPI_COMM_WORLD, atoms, gaussian_period, 2, false);
    const std::set<Key> local_gaussian_old = make_local_keys(independently_distributed_gaussian);

    unsigned long long local_common_overlap
        = static_cast<unsigned long long>(count_intersection(local_bare, local_gaussian_on_bare_owner));
    unsigned long long local_old_overlap
        = static_cast<unsigned long long>(count_intersection(local_bare, local_gaussian_old));
    unsigned long long local_bare_count = static_cast<unsigned long long>(local_bare.size());
    unsigned long long local_old_count = static_cast<unsigned long long>(local_gaussian_old.size());
    unsigned long long common_overlap = 0;
    unsigned long long old_overlap = 0;
    unsigned long long bare_count = 0;
    unsigned long long old_count = 0;
    MPI_Allreduce(&local_common_overlap, &common_overlap, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_old_overlap, &old_overlap, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_bare_count, &bare_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_old_count, &old_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    const auto expected_bare = static_cast<unsigned long long>(make_full_keys(atoms, bare_period).size());
    const auto expected_gaussian = static_cast<unsigned long long>(gaussian_support.size());
    const bool passed = bare_count == expected_bare && old_count == expected_gaussian
                        && common_overlap == expected_gaussian && old_overlap < expected_gaussian;

    if (mpi_rank == 0)
    {
        std::cout << "bare tasks: " << bare_count << "/" << expected_bare << '\n'
                  << "Gaussian tasks: " << old_count << "/" << expected_gaussian << '\n'
                  << "common-owner overlap: " << common_overlap << "/" << expected_gaussian << '\n'
                  << "old independent-owner overlap: " << old_overlap << "/" << expected_gaussian << '\n';
        if (!passed)
        {
            std::cerr << "Ewald MPI real-space owner regression failed\n";
        }
    }

    MPI_Finalize();
    return passed ? 0 : 1;
}
