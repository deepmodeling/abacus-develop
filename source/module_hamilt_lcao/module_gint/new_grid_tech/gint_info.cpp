#include <cmath>
#include <map>
#include "module_parameter/parameter.h"
#include "gint_info.h"
#include "gint_type.h"

namespace Gint
{

GintInfo::GintInfo(
    int nbx, int nby, int nbz,
    int nmx, int nmy, int nmz,
    int startidx_bx, int startidx_by, int startidx_bz,
    int nbx_local, int nby_local, int nbz_local,
    const Numerical_Orbital* Phi,
    const UnitCell& ucell, Grid_Driver& gd)
    : ucell_(&ucell)
{
    // initialize the unitcell information
    unitcell_info_ = std::make_shared<UnitCellInfo>(ucell_->a1, ucell_->a2, ucell_->a3,
                                                    nbx, nby, nbz, nmx, nmy, nmz);

    // initialize the divide information
    divide_info_ = std::make_shared<DivideInfo>(startidx_bx, startidx_by, startidx_bz,
                                                nbx, nby, nbz, unitcell_info_, false);

    // initialize the localcell information
    localcell_info_ = divide_info_->get_localcell_info();

    // initialize the biggrids
    for (int i = 0; i < localcell_info_->get_biggrid_num(); i++)
    {
        biggrids_.push_back(
            std::make_shared<BigGrid>(localcell_info_->get_bgrid_global_idx_1D(i), localcell_info_));
    }

    is_atom_in_proc_ = std::vector<bool>(ucell.nat, false);

    // initialize the atoms
    init_atoms_(ucell_->ntype, ucell_->atoms, Phi);

    // initialize the ijr_info
    init_ijr_info_(ucell, gd);
}

template <typename T>
std::shared_ptr<hamilt::HContainer<T>> GintInfo::get_hr(int npol) const
{
    auto p_hr = std::make_shared<hamilt::HContainer<T>>(ucell_->nat);
    if(PARAM.inp.gamma_only)
    {
        p_hr->fix_gamma();
    }
    p_hr->insert_ijrs(&ijr_info_, *ucell_, npol);
    p_hr->allocate(nullptr, true);
    return p_hr;
}

void GintInfo::init_atoms_(int ntype, const Atom* atoms, const Numerical_Orbital* Phi)
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

// TODO: USE OPENMP TO PARALLELIZE THIS LOOP
    for(int i = 0; i < ntype; i++)
    {
        const auto atom = atoms[i];
        const auto *orb = &Phi[i];

        // rcut extends to the maximum big grids in x, y, z directions
        int ext_bgrid_x = static_cast<int>(atom.Rcut * g1) + 1;
        int ext_bgrid_y = static_cast<int>(atom.Rcut * g2) + 1;
        int ext_bgrid_z = static_cast<int>(atom.Rcut * g3) + 1;

        for(int j = 0; j < atom.na; j++)
        {
            const Vec3i atom_bgrid_idx = unitcell_info_->get_biggrid_idx_3d(atom.tau[j]);
            const Vec3d tau_in_biggrid = atom.tau[j] - unitcell_info_->get_biggrid_coord(atom_bgrid_idx);
            std::map<Vec3i, std::shared_ptr<GintAtom>> gint_atom_map;
            for(int bgrid_x = atom_bgrid_idx.x - ext_bgrid_x; bgrid_x <= atom_bgrid_idx.x + ext_bgrid_x; bgrid_x++)
            {
                for(int bgrid_y = atom_bgrid_idx.y - ext_bgrid_y; bgrid_y <= atom_bgrid_idx.y + ext_bgrid_y; bgrid_y++)
                {
                    for(int bgrid_z = atom_bgrid_idx.z - ext_bgrid_z; bgrid_z <= atom_bgrid_idx.z + ext_bgrid_z; bgrid_z++)
                    {
                        // get the extended biggrid idx of the affected biggrid
                        const Vec3i ext_bgrid_idx(bgrid_x, bgrid_y, bgrid_z);
                        const Vec3i normal_bgrid_idx = unitcell_info_->map_ext_idx_to_ucell(ext_bgrid_idx);
                        if(localcell_info_->is_bgrid_in_lcell(normal_bgrid_idx) == false)
                        {
                            continue;
                        }
                        const int bgrid_local_idx = localcell_info_->get_bgrid_local_idx_1D(normal_bgrid_idx);
                        // get the unitcell idx of the extended atom
                        const Vec3i unitcell_idx = - unitcell_info_->get_unitcell_idx(ext_bgrid_idx);
                        auto it = gint_atom_map.find(unitcell_idx);
                        // if the gint_atom is not in the map,
                        // it means this is the first time we find this atom may affect some biggrids.
                        if(it == gint_atom_map.end())
                        {
                            Vec3i ext_atom_bgrid_idx(atom_bgrid_idx.x + unitcell_idx.x * unitcell_info_->get_nbx(),
                                                     atom_bgrid_idx.y + unitcell_idx.y * unitcell_info_->get_nby(),
                                                     atom_bgrid_idx.z + unitcell_idx.z * unitcell_info_->get_nbz()); 
                            auto gint_atom = std::make_shared<GintAtom>(&atom, j, iat, ext_atom_bgrid_idx, unitcell_idx, tau_in_biggrid, orb);
                            gint_atom_map[unitcell_idx] = gint_atom;
                            if(biggrids_[bgrid_local_idx]->is_atom_on_bgrid(*gint_atom))
                            {
                                #pragma omp critical
                                {
                                    biggrids_[bgrid_local_idx]->add_atom(gint_atom);
                                    atoms_.push_back(gint_atom);
                                    is_atom_in_proc_[iat] = true;
                                }
                            }
                        } else
                        {
                            if(biggrids_[bgrid_local_idx]->is_atom_on_bgrid(*gint_atom_map[unitcell_idx]))
                            {
                                #pragma omp critical
                                {
                                    biggrids_[bgrid_local_idx]->add_atom(gint_atom_map[unitcell_idx]);
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

void GintInfo::init_ijr_info_(const UnitCell& ucell, Grid_Driver& gd)
{
    hamilt::HContainer<double> hRGint_local(ucell.nat);

    // prepare the row_index and col_index for construct AtomPairs, they are
    // same, name as orb_index
    std::vector<int> orb_index(ucell.nat + 1);
    orb_index[0] = 0;
    for (int i = 1; i < orb_index.size(); i++) {
        int type = ucell.iat2it[i - 1];
        orb_index[i] = orb_index[i - 1] + ucell.atoms[type].nw;
    }

    for (int T1 = 0; T1 < ucell.ntype; ++T1) {
            const Atom* atom1 = &(ucell.atoms[T1]);
            for (int I1 = 0; I1 < atom1->na; ++I1) {
                auto& tau1 = atom1->tau[I1];
                const int iat1 = ucell.itia2iat(T1, I1);
                // whether this atom is in this processor.
                if (this->is_atom_in_proc_[iat1]) {
                    gd.Find_atom(ucell, tau1, T1, I1);
                    for (int ad = 0; ad < gd.getAdjacentNum() + 1; ++ad) {
                        const int T2 = gd.getType(ad);
                        const int I2 = gd.getNatom(ad);
                        const int iat2 = ucell.itia2iat(T2, I2);
                        const Atom* atom2 = &(ucell.atoms[T2]);

                        // NOTE: hRGint wil save total number of atom pairs,
                        // if only upper triangle is saved, the lower triangle will
                        // be lost in 2D-block parallelization. if the adjacent atom
                        // is in this processor.
                        if (this->is_atom_in_proc_[iat2]) {
                            Vec3d dtau = gd.getAdjacentTau(ad) - tau1;
                            double distance = dtau.norm() * ucell.lat0;
                            double rcut = atom1->Rcut + atom2->Rcut;

                            // if(distance < rcut)
                            //  mohan reset this 2013-07-02 in Princeton
                            //  we should make absolutely sure that the distance is
                            //  smaller than rcuts[it] this should be consistant
                            //  with LCAO_nnr::cal_nnrg function typical example : 7
                            //  Bohr cutoff Si orbital in 14 Bohr length of cell.
                            //  distance = 7.0000000000000000
                            //  rcuts[it] = 7.0000000000000008
                            if (distance < rcut - 1.0e-15) {
                                // calculate R index
                                auto& R_index = gd.getBox(ad);
                                // insert this atom-pair into this->hRGint
                                hamilt::AtomPair<double> tmp_atom_pair(
                                    iat1,
                                    iat2,
                                    R_index.x,
                                    R_index.y,
                                    R_index.z,
                                    orb_index.data(),
                                    orb_index.data(),
                                    ucell.nat);
                                hRGint_local.insert_pair(tmp_atom_pair);
                            }
                        }
                    }
                }
            }
    }
    this->ijr_info_ = hRGint_local.get_ijr_info();
    return;
}

template std::shared_ptr<hamilt::HContainer<double>> GintInfo::get_hr<double>(int npol) const;
template std::shared_ptr<hamilt::HContainer<std::complex<double>>> GintInfo::get_hr<std::complex<double>>(int npol) const;
}