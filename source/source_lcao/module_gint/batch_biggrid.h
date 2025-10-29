#pragma once
#include <memory>
#include <vector>
#include "big_grid.h"

namespace ModuleGint
{

class BatchBigGrid
{
    public:
    BatchBigGrid(std::vector<std::shared_ptr<BigGrid>> biggrids);
    
    const std::vector<std::shared_ptr<BigGrid>>& get_bgrids() { return biggrids_; }

    int get_batch_size() const { return biggrids_.size(); }
    int get_atoms_num() const { return atoms_num_; }
    int get_phi_len() const { return phi_len_;}
    int get_max_atoms_num_per_bgrid() const { return max_atoms_num_per_bgrid_; }
    bool empty() {return atoms_num_ == 0; }
    static int get_max_batch_size() { return max_batch_size_; }
    static int get_max_atoms_num() { return max_atoms_num_; }
    static int get_max_phi_len() { return max_phi_len_; }
    static int get_max_atom_pairs_num() { return max_atom_pairs_num_; }
    static std::shared_ptr<const BigGridInfo> get_bgrid_info() { return BigGrid::get_bgrid_info(); }
    
    private:
    std::vector<std::shared_ptr<BigGrid>> biggrids_;

    // the max nw of an atom
    int max_nw_ = 0;
    
    int phi_len_ = 0;
    // number of atoms in the batch
    int atoms_num_ = 0;

    // the max number of atoms of a single biggrid
    int max_atoms_num_per_bgrid_ = 0;

    // the max number of biggrids of a biggrids batch
    static int max_batch_size_;
    // the max number of total atoms of a biggrids batch
    static int max_atoms_num_;
    // the max number of total wavefunctions of a biggrids batch
    static int max_phi_len_;
    // the max number of atom pairs of a biggrids batch
    static int max_atom_pairs_num_;
};

}