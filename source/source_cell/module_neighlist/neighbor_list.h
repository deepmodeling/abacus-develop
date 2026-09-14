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

    void initialize(std::size_t ncentral_atoms, std::size_t pgsize)
    {
        ncentral_atoms_ = ModuleNeighList::checked_int_size(ncentral_atoms, "NeighborList central atom count");
        allocator_.initialize(ModuleNeighList::checked_int_size(pgsize, "NeighborList page size"));
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
    PageAllocator& get_allocator() { return allocator_; }
    const PageAllocator& get_allocator() const { return allocator_; }

    void set_neighbors(int i, const std::vector<int>& neighbors)
    {
        numneigh_[i] = ModuleNeighList::checked_int_size(neighbors.size(), "NeighborList neighbor count");
        firstneigh_[i] = allocator_.allocate(numneigh_[i]);
        for (int j = 0; j < numneigh_[i]; ++j)
        {
            firstneigh_[i][j] = neighbors[static_cast<std::size_t>(j)];
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
