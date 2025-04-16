#pragma once
#include <memory>
#include <vector>
#include "big_grid.h"

namespace ModuleGint
{

class BatchBigGrid
{
    public:
    BatchBigGrid(std::vector<std::shared_ptr<BigGrid>> biggrids)
    {
        biggrids_ = biggrids;
        max_batch_size_ = std::max(max_batch_size_, biggrids_.size());
        int atoms_num = 0;
        int phi_len = 0;
        for(const auto& biggrid : biggrids_)
        {
            atoms_num += biggrid->get_atoms_num();
            phi_len += biggrid->get_phi_len();
        }
        max_atoms_num_ = std::max(max_atoms_num_, atoms_num);
        max_phi_len_ = std::max(max_phi_len_, phi_len);
    };
    
    const std::vector<std::shared_ptr<BigGrid>>& get_bgrids() { return biggrids_; } 

    static int get_max_batch_size() { return max_batch_size_; }
    static int get_max_atoms_num() { return max_atoms_num_; }
    static int get_max_phi_len() { return max_phi_len_; }
    
    private:
    std::vector<std::shared_ptr<BigGrid>> biggrids_;

    // the max number of biggrids of a batch_biggrid
    static int max_batch_size_;
    // the max number of total atoms of a batch_biggrid
    static int max_atoms_num_;
    // the max number of total wavefunctions of a batch_biggrid
    static int max_phi_len_;
}

}