#include "batch_biggrid.h"

namespace ModuleGint
{

int BatchBigGrid::max_batch_size_ = 0;
int BatchBigGrid::max_atoms_num_ = 0;
int BatchBigGrid::max_phi_len_ = 0;
int BatchBigGrid::max_atom_pairs_num_ = 0;

BatchBigGrid::BatchBigGrid(std::vector<std::shared_ptr<BigGrid>> biggrids)
{
    biggrids_ = biggrids;
    max_batch_size_ = std::max(max_batch_size_, (int)biggrids_.size());
    int atom_pairs_num = 0;
    for(const auto& biggrid : biggrids_)
    {
        for(const auto& atom: biggrid->get_atoms())
        {
            max_nw_ = std::max(max_nw_, atom->get_nw());
        }
        max_atoms_num_per_bgrid_ = std::max(max_atoms_num_per_bgrid_, biggrid->get_atoms_num());
        atoms_num_ += biggrid->get_atoms_num();
        atom_pairs_num += std::pow(biggrid->get_atoms_num(), 2);
        phi_len_ += biggrid->get_phi_len() * biggrid->get_mgrids_num();
    }
    max_atoms_num_ = std::max(max_atoms_num_, atoms_num_);
    max_phi_len_ = std::max(max_phi_len_, phi_len_);
    max_atom_pairs_num_ = std::max(max_atom_pairs_num_, atom_pairs_num);
}



}
