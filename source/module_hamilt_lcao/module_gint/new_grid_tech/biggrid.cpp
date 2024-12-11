#include "biggrid.h"

namespace Gint
{

BigGrid::BigGrid(int idx, std::shared_ptr<const LocalCellInfo> localcell_info)
    : idx_(idx), localcell_info_(localcell_info),
      unitcell_info_(localcell_info->get_unitcell_info()),
      biggrid_info_(localcell_info_->get_biggrid_info()) {}

void BigGrid::add_atom(std::shared_ptr<GintAtom> atom)
{
    atoms_.push_back(atom);
}

int BigGrid::get_mgrid_psi_len() const
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

std::vector<int> BigGrid::get_atom_psi_len() const
{
    std::vector<int> psi_len(atoms_.size());
    for(int i = 0; i < atoms_.size(); ++i)
    {
        psi_len[i] = atoms_[i]->get_nw();
    }
    return psi_len;
}

std::vector<Vec3d> BigGrid::get_mgrid_coords(Vec3i bgrid_idx, Vec3d tau_in_bgrid) const
{
    Vec3i this_bgrid_idx = localcell_info_->get_bgrid_global_idx_3D(idx_);
    
    // the relative coordinates of this big grid and the atom
    Vec3d bgrid_relative_coords 
        = unitcell_info_->get_relative_coord(bgrid_idx, this_bgrid_idx) + tau_in_bgrid;

    std::vector<Vec3d> coords(biggrid_info_->get_nmxyz());
    for(int im = 0; im < biggrid_info_->get_nmxyz(); ++im)
    {
        Vec3d mcell_coords = biggrid_info_->get_meshgrid_coord(im);
        coords[im] = mcell_coords - bgrid_relative_coords;
    }
    return coords;
}


std::vector<Vec3d> BigGrid::get_mgrid_coords(const GintAtom& atom) const
{
    return get_mgrid_coords(atom.get_biggrid_idx(), atom.get_tau_in_biggrid());
}

std::vector<std::vector<bool>> BigGrid::get_is_atom_on_mgrid() const
{
    std::vector<std::vector<bool>> if_on_mcell;
    for(const auto& atom : atoms_)
    {
        std::vector<Vec3d> mcell_coords = get_mgrid_coords(*atom);
        std::vector<bool> atom_if_on_mcell;
        for(const auto& coord : mcell_coords)
        {
            coord.norm() < atom->get_rcut() ? atom_if_on_mcell.push_back(true) : atom_if_on_mcell.push_back(false);
        }
        if_on_mcell.push_back(atom_if_on_mcell);
    }
    return ;
}

void BigGrid::get_psi(const std::vector<Vec3d> coords, double * psi) const
{
    int psi_len = get_mgrid_psi_len();
    for(const auto& atom : atoms_)
    {
        atom->get_psi(coords, psi_len, psi);
        psi += atom->get_nw();
    }
}

void BigGrid::get_psi_dpsir(
    const std::vector<Vec3d> coords,
    double* psi, double* dpsi_x, double* dpsi_y, double* dpsi_z) const
{
    int psi_len = get_mgrid_psi_len();
    for(const auto& atom : atoms_)
    {
        atom->get_psi_dpsir(coords, psi_len, psi, dpsi_x, dpsi_y, dpsi_z);
        psi += atom->get_nw();
        dpsi_x += atom->get_nw();
        dpsi_y += atom->get_nw();
        dpsi_z += atom->get_nw();
    }
}

void BigGrid::get_ddpsir(
    const std::vector<Vec3d> coords,
    double* ddpsi_xx, double* ddpsi_xy, double* ddpsi_xz,
    double* ddpsi_yy, double* ddpsi_yz, double* ddpsi_zz) const
{
    int psi_len = get_mgrid_psi_len();
    for(const auto& atom : atoms_)
    {
        atom->get_ddpsir(coords, psi_len, ddpsi_xx, ddpsi_xy, ddpsi_xz, ddpsi_yy, ddpsi_yz, ddpsi_zz);
        ddpsi_xx += atom->get_nw();
        ddpsi_xy += atom->get_nw();
        ddpsi_xz += atom->get_nw();
        ddpsi_yy += atom->get_nw();
        ddpsi_yz += atom->get_nw();
        ddpsi_zz += atom->get_nw();
    }
}


} // namespace Gint