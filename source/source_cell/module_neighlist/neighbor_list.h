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

    void initialize(std::size_t nlocal, std::size_t pgsize)
    {
        nlocal_ = ModuleNeighList::checked_int_size(nlocal, "NeighborList local atom count");
        allocator_ = PageAllocator(ModuleNeighList::checked_int_size(pgsize, "NeighborList page size"));
        numneigh_.assign(nlocal, 0);
        firstneigh_.assign(nlocal, nullptr);
    }

    void reset()
    {
        allocator_.reset();
    }

    int get_nlocal() const { return nlocal_; }
    int get_numneigh(int i) const { return numneigh_[i]; }
    int* get_firstneigh(int i) { return firstneigh_[i]; }
    const int* get_firstneigh(int i) const { return firstneigh_[i]; }
    PageAllocator& get_allocator() { return allocator_; }
    const PageAllocator& get_allocator() const { return allocator_; }

private:
    int nlocal_ = 0;
    std::vector<int> numneigh_;
    std::vector<int*> firstneigh_;
    PageAllocator allocator_;

    friend class BinManager;
};

#endif // NEIGHBOR_LIST_H
