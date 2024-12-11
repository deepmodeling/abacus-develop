#include <cmath>
#include <unordered_map>
#include "gint_info.h"
#include "gint_type.h"

namespace Gint
{

GintInfo::GintInfo(
    Vec3d unitcell_vec1,
    Vec3d unitcell_vec2,
    Vec3d unitcell_vec3,
    int nbx, int nby, int nbz,
    int nmx, int nmy, int nmz,
    int startidx_bx, int startidx_by, int startidx_bz,
    int nbx_local, int nby_local, int nbz_local,
    int ntype, Atom* atoms, Numerical_Orbital* Phi)
{
    // initialize the unitcell information
    unitcell_info_ = std::make_shared<UnitCellInfo>(unitcell_vec1, unitcell_vec2, unitcell_vec3,
                                                    nbx, nby, nbz, nmx, nmy, nmz);

    // initialize the divide information
    divide_info_ = std::make_shared<DivideInfo>(startidx_bx, startidx_by, startidx_bz,
                                                nbx, nby, nbz, unitcell_info_, false);

    // initialize the localcell information
    localcell_info_ = divide_info_->get_localcell_info();

    // initialize the biggrids
    for (int i = 0; i < localcell_info_->get_biggrid_num(); i++)
    {
        biggrids_.emplace_back(
            std::make_shared<BigGrid>(localcell_info_->get_bgrid_global_idx_1D(i), localcell_info_));
    }

    // initialize the atoms
    init_atoms_(ntype, atoms, Phi);
}

void GintInfo::init_atoms_(int ntype, Atom* atoms, Numerical_Orbital* Phi)
{
    int iat = 0;
    const Matrix3 biggrid_GT = unitcell_info_->get_biggrid_info()->get_GT();
    const double g1 = sqrt(biggrid_GT.e11 * biggrid_GT.e11
        + biggrid_GT.e21 * biggrid_GT.e21
        + biggrid_GT.e31 * biggrid_GT.e31);
    const double g2 = sqrt(biggrid_GT.e12 * biggrid_GT.e12
        + biggrid_GT.e22 * biggrid_GT.e22
        + biggrid_GT.e32 * biggrid_GT.e32);
    const double g3 = sqrt(biggrid_GT.e13 * biggrid_GT.e13
        + biggrid_GT.e23 * biggrid_GT.e23
        + biggrid_GT.e33 * biggrid_GT.e33);

#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < ntype; i++)
    {
        const auto atom = atoms[i];
        const auto *orb = &Phi[i];
        // rcut extends to the maximum big grids in x, y, z directions
        int rcut_bgrid_x = static_cast<int>(atom.Rcut * g1) + 1;
        int rcut_bgrid_y = static_cast<int>(atom.Rcut * g2) + 1;
        int rcut_bgrid_z = static_cast<int>(atom.Rcut * g3) + 1;

        double rcut_ucell_x = rcut_bgrid_x / (double)unitcell_info_->get_nbx();
        double rcut_ucell_y = rcut_bgrid_y / (double)unitcell_info_->get_nby();
        double rcut_ucell_z = rcut_bgrid_z / (double)unitcell_info_->get_nbz();

        for(int j = 0; j < atom.na; j++)
        {
            const Vec3i biggrid_idx = unitcell_info_->get_biggrid_idx_3d(atom.tau[j]);
            const Vec3d tau_in_biggrid = atom.tau[j] - unitcell_info_->get_biggrid_coord(biggrid_idx);
            int min_x = static_cast<int>(std::floor(atom.taud[j].x - rcut_ucell_x));
            int max_x = static_cast<int>(std::floor(atom.taud[j].x + rcut_ucell_x));
            int min_y = static_cast<int>(std::floor(atom.taud[j].y - rcut_ucell_y));
            int max_y = static_cast<int>(std::floor(atom.taud[j].y + rcut_ucell_y));
            int min_z = static_cast<int>(std::floor(atom.taud[j].z - rcut_ucell_z));
            int max_z = static_cast<int>(std::floor(atom.taud[j].z + rcut_ucell_z));
            for(int k = min_x; k <= max_x; k++)
            {
                for(int l = min_y; l <= max_y; l++)
                {
                    for(int m = min_z; m <= max_z; m++)
                    {
                        const Vec3i atom_idx_ext(biggrid_idx.x + k * unitcell_info_->get_nbx(),
                                                biggrid_idx.y + l * unitcell_info_->get_nby(),
                                                biggrid_idx.z + m * unitcell_info_->get_nbz());
                        const Vec3i unitcell_idx(k, l, m);
                        auto gint_atom = std::make_shared<GintAtom>(&atom, iat, atom_idx_ext, unitcell_idx, tau_in_biggrid, orb);
                        #pragma omp critical
                        { atoms_.push_back(gint_atom); }
                        for(int bx = atom_idx_ext.x - rcut_bgrid_x; bx <= atom_idx_ext.x + rcut_bgrid_x; bx++)
                        {
                            for(int by = atom_idx_ext.y - rcut_bgrid_y; by <= atom_idx_ext.y + rcut_bgrid_y; by++)
                            {
                                for(int bz = atom_idx_ext.z - rcut_bgrid_z; bz <= atom_idx_ext.z + rcut_bgrid_z; bz++)
                                {
                                    const Vec3i biggrid_idx = unitcell_info_->map_ext_idx_to_ucell(Vec3i(bx, by, bz));
                                    if (localcell_info_->is_bgrid_in_lcell(biggrid_idx))
                                    {   
                                        const auto biggrid = biggrids_[localcell_info_->get_bgrid_local_idx_1D(biggrid_idx)];
                                        const auto meshgrid_coords = biggrid->get_mgrid_coords(*gint_atom);
                                        for(const auto& coord : meshgrid_coords)
                                        {
                                            if(coord.norm() < atom.Rcut)
                                            {
                                                #pragma omp critical
                                                { biggrid->add_atom(gint_atom);}
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            iat++;
        }
    }
}
}