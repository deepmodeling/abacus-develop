#pragma once

#include <vector>
#include <memory>
#include "gint_type.h"
#include "biggrid_info.h"
#include "localcell_info.h"
#include "unitcell_info.h"
#include "gint_atom.h"

namespace ModuleGint
{

class BigGrid
{
    public:
        // constructor
        BigGrid(int idx);

        static void init_localcell_info(std::shared_ptr<const LocalCellInfo> localcell_info) { localcell_info_ = localcell_info; };
        static void init_unitcell_info(std::shared_ptr<const UnitCellInfo> unitcell_info) { unitcell_info_ = unitcell_info; };
        static void init_bgrid_info(std::shared_ptr<const BigGridInfo> biggrid_info) { biggrid_info_ = biggrid_info; };

        // getter functions
        int get_idx() const { return idx_; };
        std::shared_ptr<const LocalCellInfo> get_localcell_info() const { return localcell_info_; };
        std::shared_ptr<const UnitCellInfo> get_unitcell_info() const {return unitcell_info_; };
        std::shared_ptr<const BigGridInfo> get_bgrid_info() const { return biggrid_info_; };
        const std::vector<const GintAtom*>& get_atoms() const { return atoms_; };
        const GintAtom* get_atom(int i) const { return atoms_[i]; };

        // get the number of meshgrids in the big grid
        int get_mgrids_num() const { return biggrid_info_->get_mgrids_num(); };

        // get the number of atoms that can affect the big grid
        int get_atoms_num() const { return atoms_.size(); };

        // add an atom to the big grid
        void add_atom(const GintAtom* atom);

        // get the total number of phi of a meshgrid
        // return: (\sum_{i=0}^{atoms_->size()} atoms_[i]->nw)
        int get_mgrid_phi_len() const;

        // set the start index of the phi of each atom
        // return: vector[i] = \sum_{j=0}^{i-1} atoms_[j]->nw
        void set_atoms_startidx(std::vector<int>& startidx) const;

        // set the length of phi of each atom
        void set_atoms_phi_len(std::vector<int>& phi_len) const;

        // set the coordinates of the meshgrids of the big grid
        void set_mgrids_coord(std::vector<Vec3d>& coord) const;

        // set the 1D index of the meshgrids in the local cell
        void set_mgrids_local_idx(std::vector<int>& mgrids_idx) const;

        /**
         * @brief Set the coordinates of the meshgrids of the big grid relative to an atom
         * 
         * @param bgrid_idx the 3D index of the big grid, which contains the atom, in the unitcell
         * @param tau_in_bgrid the cartesian coordinate of the atom relative to the big grid containing it
         * @param atom_coord the relative cartesian coordinates of the atom and the meshgrids
         */
        void set_atom_relative_coords(const Vec3i bgrid_idx, const Vec3d tau_in_bgrid, std::vector<Vec3d>& atom_coord) const;

        // a wrapper function to get the relative coordinates of the atom and the meshgrids
        void set_atom_relative_coords(const GintAtom* atom, std::vector<Vec3d>& atom_coord) const;

        // if the atom affects the big grid, return true, otherwise false
        // note when we say an atom affects a big grid, it does not mean that the atom affects all the meshgrid on the big grid,
        // it may only affect a part of them.
        bool is_atom_on_bgrid(const GintAtom* atom) const;
    
    private:
        // atoms that can affect the big grid
        std::vector<const GintAtom*> atoms_;

        // the 1D index of the big grid in the local cell
        const int idx_;

        // local cell info
        static std::shared_ptr<const LocalCellInfo> localcell_info_;

        // unitcell info
        static std::shared_ptr<const UnitCellInfo> unitcell_info_;

        // the big grid info
        static std::shared_ptr<const BigGridInfo> biggrid_info_;
};

} // namespace ModuleGint