#pragma once

#include <vector>
#include <memory>
#include "gint_type.h"
#include "biggrid_info.h"
#include "localcell_info.h"
#include "unitcell_info.h"
#include "gint_atom.h"

namespace Gint
{

class BigGrid
{
    public:
        // constructor
        BigGrid(int idx, std::shared_ptr<const LocalCellInfo> localcell_info);

        // getter functions
        int get_idx() const { return idx_; };
        std::shared_ptr<const LocalCellInfo> get_localcell_info() const { return localcell_info_; };
        std::shared_ptr<const UnitCellInfo> get_unitcell_info() const {return unitcell_info_; };
        std::shared_ptr<const BigGridInfo> get_biggrid_info() const { return biggrid_info_; };
        const std::vector<std::shared_ptr<GintAtom>>& get_atoms() const { return atoms_; };

        // get the number of meshgrids in the big grid
        int get_meshgrid_num() const { return biggrid_info_->get_nmxyz(); };

        // get the number of atoms that can affect the big grid
        int get_atom_num() const { return atoms_.size(); };

        // add an atom to the big grid
        void add_atom(std::shared_ptr<GintAtom> atom);

        // get the total number of phi of a meshgrid
        // return: (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
        int get_mgrid_phi_len() const;

        // get the start index of the phi of each atom
        // return: vector[i] = \sum_{j=0}^{i-1} atoms_[j]->nw
        vector<int> get_atom_startidx() const;

        // get the length of phi of each atom
        vector<int> get_atom_phi_len() const;

        // get the coordinates of the meshgrids of the big grid
        vector<Vec3d> get_mgrid_coords() const;

        /**
         * @brief Get the coordinates of the meshgrids of the big grid relative to an atom
         * 
         * @param bgrid_idx the 3D index of the big grid, which contains the atom, in the unitcell
         * @param tau_in_bgrid the cartesian coordinate of the atom relative to the big grid containing it
         */
        std::vector<Vec3d> get_atom_relative_coords(Vec3i bgrid_idx, Vec3d tau_in_bgrid) const;

        // get the coordinates of the meshgrids of the big grid relative to the atom
        std::vector<Vec3d> get_atom_relative_coords(const GintAtom& atom) const;

        // if the atom affects the big grid, return true, otherwise false
        // note when we say an atom affects a big grid, it does not mean that the atom affects all the meshgrid on the big grid,
        // it may only affect a part of them.
        bool is_atom_on_bgrid(const GintAtom& atom) const;

        /**
         * @brief Get a boolean array to indicate whether the atom affects the meshgrids
         * 
         * @return a boolean array, the dimension of the array is atoms_->size() * biggrid_info_->get_nmxyz()
         * array[i][j] = true if the ith atom affects the jth meshgrid，otherwise false.
         */
        std::vector<std::vector<bool>> get_is_atom_on_mgrids() const;
    
    private:
        // atoms that can affect the big grid
        std::vector<std::shared_ptr<GintAtom>> atoms_;

        // the 1D index of the big grid in the local cell
        const int idx_;

        // local cell info
        std::shared_ptr<const LocalCellInfo> localcell_info_;

        // unitcell info
        std::shared_ptr<const UnitCellInfo> unitcell_info_;

        // the big grid info
        std::shared_ptr<const BigGridInfo> biggrid_info_;
};

} // namespace Gint