#include "big_grid.h"

namespace ModuleGint
{
std::shared_ptr<const LocalCellInfo> BigGrid::localcell_info_ = nullptr;
std::shared_ptr<const UnitCellInfo> BigGrid::unitcell_info_ = nullptr;
std::shared_ptr<const BigGridInfo> BigGrid::biggrid_info_ = nullptr;

BigGrid::BigGrid(int idx): idx_(idx){}

void BigGrid::add_atom(const GintAtom* atom)
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

void BigGrid::set_atoms_startidx(std::vector<int>& startidx) const
{
    startidx.resize(atoms_.size());
    startidx[0] = 0;
    for(int i = 1; i < atoms_.size(); ++i)
    {
        startidx[i] = startidx[i-1] + atoms_[i-1]->get_nw();
    }
}

void BigGrid::set_atoms_phi_len(std::vector<int>& phi_len) const
{
    phi_len.resize(atoms_.size());
    for(int i = 0; i < atoms_.size(); ++i)
    {
        phi_len[i] = atoms_[i]->get_nw();
    }
}

void BigGrid::set_mgrids_coord(std::vector<Vec3d>& coord) const
{
    coord.resize(biggrid_info_->get_mgrids_num());
    Vec3d this_bgrid_coord = localcell_info_->get_bgrid_global_coord_3D(idx_);
    for(int im = 0; im < biggrid_info_->get_mgrids_num(); ++im)
    {
        coord[im] = biggrid_info_->get_mgrid_coord(im) + this_bgrid_coord;
    }
}

void BigGrid::set_mgrids_local_idx(std::vector<int>& mgrids_idx) const
{
    auto index_3d = localcell_info_->bgrid_idx_1Dto3D(idx_);
    Vec3i startidx(
        index_3d.x * biggrid_info_->get_nmx(),
        index_3d.y * biggrid_info_->get_nmy(),
        index_3d.z * biggrid_info_->get_nmz());
    mgrids_idx.resize(0);
    for(int ix = 0; ix < biggrid_info_->get_nmx(); ++ix)
    {
        for(int iy = 0; iy < biggrid_info_->get_nmy(); ++iy)
        {
            for(int iz = 0; iz < biggrid_info_->get_nmz(); ++iz)
            {
                Vec3i idx_3d(startidx.x + ix, startidx.y + iy, startidx.z + iz);
                mgrids_idx.push_back(localcell_info_->mgrid_idx_3Dto1D(idx_3d));
            }
        }
    }
}

void BigGrid::set_atom_relative_coords(const Vec3i bgrid_idx, const Vec3d tau_in_bgrid, std::vector<Vec3d>& atom_coord) const
{
    Vec3i this_bgrid_idx = localcell_info_->get_bgrid_global_idx_3D(idx_);
    
    // the relative coordinates of this big grid and the atom
    Vec3d bgrid_relative_coord 
        = unitcell_info_->get_relative_coord(bgrid_idx, this_bgrid_idx) + tau_in_bgrid;

    atom_coord.resize(biggrid_info_->get_mgrids_num());
    for(int im = 0; im < biggrid_info_->get_mgrids_num(); ++im)
    {
        const Vec3d& mcell_coord = biggrid_info_->get_mgrid_coord(im);
        atom_coord[im] = mcell_coord - bgrid_relative_coord;
    }
}


void BigGrid::set_atom_relative_coords(const GintAtom* atom, std::vector<Vec3d>& atom_coord) const
{
    set_atom_relative_coords(atom->get_bgrid_idx(), atom->get_tau_in_bgrid(), atom_coord);
}

bool BigGrid::is_atom_on_bgrid(const GintAtom* atom) const
{
    std::vector<Vec3d> coords;
    this->set_atom_relative_coords(atom, coords);
    for(const auto& dist : coords)
    {
        if(dist.norm() <= atom->get_rcut())
        {
            return true;
        }
    }
    return false;
}

} // namespace ModuleGint