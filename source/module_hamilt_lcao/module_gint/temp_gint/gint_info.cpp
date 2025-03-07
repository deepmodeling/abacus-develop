#include <cmath>
#include <map>
#include "module_parameter/parameter.h"
#include "module_base/timer.h"
#include "gint_info.h"
#include "gint_type.h"

namespace ModuleGint
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
    unitcell_info_ = std::make_shared<UnitCellInfo>(ucell_->a1 * ucell_->lat0, ucell_->a2 * ucell_->lat0, ucell_->a3 * ucell_->lat0,
                                                    nbx, nby, nbz, nmx, nmy, nmz);

    biggrid_info_ = unitcell_info_->get_bgrid_info();
    meshgrid_info_ = biggrid_info_->get_mgrid_info();

    // initialize the divide information
    divide_info_ = std::make_shared<DivideInfo>(startidx_bx, startidx_by, startidx_bz,
                                                nbx_local, nby_local, nbz_local, unitcell_info_, false);

    // initialize the localcell information
    localcell_info_ = divide_info_->get_localcell_info();

    // initialize the biggrids
    BigGrid::init_localcell_info(localcell_info_);
    BigGrid::init_unitcell_info(unitcell_info_);
    BigGrid::init_bgrid_info(biggrid_info_);

    for (int i = 0; i < localcell_info_->get_bgrids_num(); i++)
    {
        biggrids_.push_back(std::make_shared<BigGrid>(i));
    }

    // initialize the atoms
    init_atoms_(ucell_->ntype, ucell_->atoms, Phi);

    // initialize the ijr_info
    // this step needs to be done after init_atoms_, because it requires the information of is_atom_on_bgrid
    init_ijr_info_(ucell, gd);
}

template <typename T>
std::shared_ptr<HContainer<T>> GintInfo::get_hr(int npol) const
{
    auto p_hr = std::make_shared<HContainer<T>>(ucell_->nat);
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
    ModuleBase::timer::tick("GintInfo", "init_atoms");
    int iat = 0;
    is_atom_in_proc_.resize(ucell_->nat, false);
    atoms_.resize(ucell_->nat);

// TODO: USE OPENMP TO PARALLELIZE THIS LOOP
    for(int i = 0; i < ntype; i++)
    {
        const auto& atom = atoms[i];
        const auto *orb = &Phi[i];

        // rcut extends to the maximum big grids in x, y, z directions
        Vec3i ext_bgrid = biggrid_info_->max_ext_bgrid_num(atom.Rcut);

        for(int j = 0; j < atom.na; j++)
        {
            Vec3d fraction;
            fraction.x = atom.taud[j].x * unitcell_info_->get_nbx();
			fraction.y = atom.taud[j].y * unitcell_info_->get_nby();
			fraction.z = atom.taud[j].z * unitcell_info_->get_nbz();
            const Vec3i atom_bgrid_idx(static_cast<int>(fraction.x),
                                       static_cast<int>(fraction.y),
                                       static_cast<int>(fraction.z));
            const Vec3d delta(fraction.x - atom_bgrid_idx.x,
                              fraction.y - atom_bgrid_idx.y,
                              fraction.z - atom_bgrid_idx.z);
            const Vec3d tau_in_biggrid = biggrid_info_->get_cartesian_coord(delta);

            const Vec3i ucell_idx_atom = unitcell_info_->get_unitcell_idx(atom_bgrid_idx);
            auto& r_to_atom = atoms_[iat];

            for(int bgrid_x = atom_bgrid_idx.x - ext_bgrid.x; bgrid_x <= atom_bgrid_idx.x + ext_bgrid.x; bgrid_x++)
            {
                for(int bgrid_y = atom_bgrid_idx.y - ext_bgrid.y; bgrid_y <= atom_bgrid_idx.y + ext_bgrid.y; bgrid_y++)
                {
                    for(int bgrid_z = atom_bgrid_idx.z - ext_bgrid.z; bgrid_z <= atom_bgrid_idx.z + ext_bgrid.z; bgrid_z++)
                    {
                        // get the extended biggrid idx of the affected biggrid
                        const Vec3i ext_bgrid_idx(bgrid_x, bgrid_y, bgrid_z);
                        const Vec3i normal_bgrid_idx = unitcell_info_->map_ext_idx_to_ucell(ext_bgrid_idx);
                        if(localcell_info_->is_bgrid_in_lcell(normal_bgrid_idx) == false)
                        {
                            continue;
                        }
                        const int bgrid_local_idx = localcell_info_->get_bgrid_local_idx_1D(normal_bgrid_idx);
                        // get the unitcell idx of the big grid
                        const Vec3i ucell_idx_bgrid = unitcell_info_->get_unitcell_idx(ext_bgrid_idx);

                        // The index of the unitcell containing the biggrid relative to the unitcell containing the atom.
                        const Vec3i ucell_idx_relative = ucell_idx_bgrid - ucell_idx_atom;
                        auto it = r_to_atom.find(ucell_idx_relative);
                        // if the gint_atom is not in the map,
                        // it means this is the first time we find this atom may affect some biggrids,
                        // add it to the r_to_atom map
                        if(it == r_to_atom.end())
                        {
                            Vec3i ext_atom_bgrid_idx(atom_bgrid_idx.x - ucell_idx_bgrid.x * unitcell_info_->get_nbx(),
                                                     atom_bgrid_idx.y - ucell_idx_bgrid.y * unitcell_info_->get_nby(),
                                                     atom_bgrid_idx.z - ucell_idx_bgrid.z * unitcell_info_->get_nbz());
                            r_to_atom.insert(std::make_pair(ucell_idx_relative, 
                                GintAtom(&atom, j, iat, ext_atom_bgrid_idx, ucell_idx_relative, tau_in_biggrid, orb)));
                        }
                        if(biggrids_[bgrid_local_idx]->is_atom_on_bgrid(&r_to_atom.at(ucell_idx_relative)))
                        {
                            biggrids_[bgrid_local_idx]->add_atom(&r_to_atom.at(ucell_idx_relative));
                            is_atom_in_proc_[iat] = true;
                        }
                    }
                }
            }
            iat++;
        }
    }
    ModuleBase::timer::tick("GintInfo", "init_atoms");
}

void GintInfo::init_ijr_info_(const UnitCell& ucell, Grid_Driver& gd)
{
    HContainer<double> hr_gint_local(ucell.nat);
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

                        // NOTE: hr_gint wil save total number of atom pairs,
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
                                // insert this atom-pair into this->hr_gint
                                hamilt::AtomPair<double> tmp_atom_pair(
                                    iat1,
                                    iat2,
                                    R_index.x,
                                    R_index.y,
                                    R_index.z,
                                    orb_index.data(),
                                    orb_index.data(),
                                    ucell.nat);
                                hr_gint_local.insert_pair(tmp_atom_pair);
                            }
                        }
                    }
                }
            }
    }
    this->ijr_info_ = hr_gint_local.get_ijr_info();
    return;
}

template std::shared_ptr<HContainer<double>> GintInfo::get_hr<double>(int npol) const;
template std::shared_ptr<HContainer<std::complex<double>>> GintInfo::get_hr<std::complex<double>>(int npol) const;
}