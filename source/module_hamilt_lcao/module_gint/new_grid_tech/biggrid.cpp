#include "biggrid.h"

namespace ModuleGint
{

BigGrid::BigGrid(int idx, std::shared_ptr<const LocalCellInfo> localcell_info)
    : idx_(idx), localcell_info_(localcell_info),
      unitcell_info_(localcell_info->get_unitcell_info()),
      biggrid_info_(localcell_info_->get_biggrid_info()) {}

void BigGrid::add_atom(std::shared_ptr<GintAtom> atom)
{
    atoms_.push_back(atom);
}

int BigGrid::get_mgrid_phi_len() const
{
    int len = 0;
    for(const auto& atom : atoms_)
    {
        len += atom->get_nw();
    }
    return len;
}

std::vector<int> BigGrid::get_atom_startidx() const
{
    std::vector<int> startidx(atoms_.size());
    startidx[0] = 0;
    for(int i = 1; i < atoms_.size(); ++i)
    {
        startidx[i] = startidx[i-1] + atoms_[i-1]->get_nw();
    }
    return startidx;
}

std::vector<int> BigGrid::get_atom_phi_len() const
{
    std::vector<int> phi_len(atoms_.size());
    for(int i = 0; i < atoms_.size(); ++i)
    {
        phi_len[i] = atoms_[i]->get_nw();
    }
    return phi_len;
}

std::vector<Vec3d> BigGrid::get_mgrid_coords() const
{
    std::vector<Vec3d> coords(biggrid_info_->get_nmxyz());
    Vec3d this_bgrid_coord = localcell_info_->get_bgrid_global_coord_3D(idx_);
    for(int im = 0; im < biggrid_info_->get_nmxyz(); ++im)
    {
        coords[im] = biggrid_info_->get_meshgrid_coord(im) + this_bgrid_coord;
    }
    return coords;
}

std::vector<int> BigGrid::get_mgrids_local_idx() const
{
    return localcell_info_->get_mgrids_local_idx_1D(idx_);
}

std::vector<Vec3d> BigGrid::get_atom_relative_coords(Vec3i bgrid_idx, Vec3d tau_in_bgrid) const
{
    Vec3i this_bgrid_idx = localcell_info_->get_bgrid_global_idx_3D(idx_);
    
    // the relative coordinates of this big grid and the atom
    Vec3d bgrid_relative_coord 
        = unitcell_info_->get_relative_coord(bgrid_idx, this_bgrid_idx) + tau_in_bgrid;

    std::vector<Vec3d> coords(biggrid_info_->get_nmxyz());
    for(int im = 0; im < biggrid_info_->get_nmxyz(); ++im)
    {
        Vec3d mcell_coord = biggrid_info_->get_meshgrid_coord(im);
        coords[im] = mcell_coord - bgrid_relative_coord;
    }
    return coords;
}


std::vector<Vec3d> BigGrid::get_atom_relative_coords(const GintAtom& atom) const
{
    return get_atom_relative_coords(atom.get_biggrid_idx(), atom.get_tau_in_biggrid());
}

bool BigGrid::is_atom_on_bgrid(const GintAtom& atom) const
{
    for(const auto& dist : get_atom_relative_coords(atom))
    {
        if(dist.norm() <= atom.get_rcut())
        {
            return true;
        }
    }
    return false;
}

std::vector<std::vector<bool>> BigGrid::get_is_atom_on_mgrids() const
{
    std::vector<std::vector<bool>> if_on_mcell;
    for(const auto& atom : atoms_)
    {
        std::vector<Vec3d> mcell_coords = get_atom_relative_coords(*atom);
        std::vector<bool> atom_if_on_mcell(mcell_coords.size());
        for(const auto& coord : mcell_coords)
        {
            for(int i = 0; i < mcell_coords.size(); ++i)
            {   
                auto coord = mcell_coords[i];
                atom_if_on_mcell[i] = coord.norm() <= atom->get_rcut();
            }
        }
        if_on_mcell.push_back(atom_if_on_mcell);
    }
    return if_on_mcell;
}

} // namespace ModuleGint