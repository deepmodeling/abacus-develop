#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H

#include "source_cell/module_neighlist/neighbor_types.h"

#include <vector>
#include "page_allocator.h"

class NeighborList
{
public:
    NeighborList() = default;
    ~NeighborList() = default;

    void initialize(int ncentral_atoms, int pgsize)
    {
        ncentral_atoms_ = ncentral_atoms;
        allocator_.initialize(pgsize);
        numneigh_.assign(ncentral_atoms, 0);
        firstneigh_.assign(ncentral_atoms, nullptr);
    }

    void reset()
    {
        allocator_.reset();
    }

    int get_ncentral_atoms() const { return ncentral_atoms_; }
    int get_numneigh(int i) const { return numneigh_[i]; }
    int* get_firstneigh(int i) { return firstneigh_[i]; }
    const int* get_firstneigh(int i) const { return firstneigh_[i]; }
    /// Return a distance-sorted copy of one atom's neighbors without changing the stored order.
    std::vector<int> get_neighbors_sorted_by_distance(int central_atom,
                                                      const std::vector<double>& coordinates) const;
    PageAllocator& get_allocator() { return allocator_; }
    const PageAllocator& get_allocator() const { return allocator_; }

    void set_neighbors(int i, const std::vector<int>& neighbors)
    {
        numneigh_[i] = static_cast<int>(neighbors.size());
        firstneigh_[i] = allocator_.allocate(numneigh_[i]);
        for (int j = 0; j < numneigh_[i]; ++j)
        {
            firstneigh_[i][j] = neighbors[j];
        }
    }

private:
    int ncentral_atoms_ = 0;
    std::vector<int> numneigh_;
    std::vector<int*> firstneigh_;
    PageAllocator allocator_;

    friend class BinManager;
};

#endif // NEIGHBOR_LIST_H
