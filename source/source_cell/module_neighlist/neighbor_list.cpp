#include "neighbor_list.h"

#include <algorithm>
#include <utility>

std::vector<int> NeighborList::get_neighbors_sorted_by_distance(const int central_atom,
                                                                const std::vector<double>& coordinates) const
{
    const int neighbor_count = get_numneigh(central_atom);
    std::vector<std::pair<double, int> > distance_neighbors;
    distance_neighbors.reserve(static_cast<std::size_t>(neighbor_count));
    for (int ineighbor = 0; ineighbor < neighbor_count; ++ineighbor)
    {
        const int neighbor_atom = get_firstneigh(central_atom)[ineighbor];
        const double dx = coordinates[3 * central_atom] - coordinates[3 * neighbor_atom];
        const double dy = coordinates[3 * central_atom + 1] - coordinates[3 * neighbor_atom + 1];
        const double dz = coordinates[3 * central_atom + 2] - coordinates[3 * neighbor_atom + 2];
        const double distance2 = dx * dx + dy * dy + dz * dz;
        distance_neighbors.push_back(std::make_pair(distance2, neighbor_atom));
    }

    std::sort(distance_neighbors.begin(), distance_neighbors.end());

    std::vector<int> sorted_neighbors;
    sorted_neighbors.reserve(static_cast<std::size_t>(neighbor_count));
    for (const std::pair<double, int>& distance_neighbor : distance_neighbors)
    {
        sorted_neighbors.push_back(distance_neighbor.second);
    }
    return sorted_neighbors;
}
