#pragma once

#include "module_cell/atom_spec.h"
#include "module_basis/module_ao/ORB_atomic.h"
#include "gint_type.h"

namespace ModuleGint
{

class GintAtom
{
    public:
        // constructor
        GintAtom(
            const Atom* atom,
            int ia,
            int iat,
            Vec3i biggrid_idx,
            Vec3i unitcell_idx,
            Vec3d tau_in_biggrid,
            const Numerical_Orbital* orb)
            : atom_(atom), ia_(ia), iat_(iat), biggrid_idx_(biggrid_idx),
              unitcell_idx_(unitcell_idx), tau_in_biggrid_(tau_in_biggrid),
              orb_(orb) {};

        // getter functions
        const Atom* get_atom() const { return atom_; };
        const int get_ia() const { return ia_; };
        const int get_iat() const { return iat_; };
        const Vec3i& get_bgrid_idx() const { return biggrid_idx_; };
        const Vec3i& get_unitcell_idx() const { return unitcell_idx_; };
        const Vec3i& get_R() const { return unitcell_idx_; };
        const Vec3d& get_tau_in_bgrid() const { return tau_in_biggrid_; };
        const Numerical_Orbital* get_orb() const { return orb_; };

        const int get_nw() const { return atom_->nw; };
        const double get_rcut() const { return orb_->getRcut(); };
        
        /**
         * @brief Get the wave function values of the atom at a meshgrid.
         * 
         * phi[(n-1)*stride] ~ phi[(n-1)*stride + nw] store the wave function values of the first atom at the nth meshgrid
         * 
         * @param coords the cartesian coordinates of the meshgrids of a biggrid relative to the atom
         * @param stride the stride of the phi array between two adjacent meshgrids
         * @param phi array to store the wave function values
         */
        template <typename T>
        void set_phi(const std::vector<Vec3d>& coords, const int stride, T* phi) const;

        /**
         * @brief Get the wave function values and its derivative
         * 
         * The reason for combining the functions to solve the wave function values 
         * and wave function derivatives into one function is to improve efficiency.
         * phi[(n-1)*stride] ~ phi[(n-1)*stride + nw] store the wave function values of the first atom at the nth meshgrid
         * 
         * @param coords the cartesian coordinates of the meshgrids of a biggrid relative to the atom
         * @param stride the stride of the phi array between two adjacent meshgrids
         * @param phi array to store the wave function values
         * @param dphi_x array to store the derivative wave functions in x direction
         * @param dphi_y array to store the derivative wave functions in y direction
         * @param dphi_z array to store the derivative wave functions in z direction
         */
        template <typename T>
        void set_phi_dphi(
            const std::vector<Vec3d>& coords, const int stride,
            T* phi, T* dphi_x, T* dphi_y, T* dphi_z) const;

        /**
         * @brief Get the wave function values and its second derivative
         * 
         * ddphi[(n-1)*stride] ~ ddphi[(n-1)*stride + nw] store the second derivative of 
         * wave function values of the atom at the first meshgrid
         *  
         * @param coords the cartesian coordinates of the meshgrids of a biggrid relative to the atom
         * @param stride the stride of the phi array between two adjacent meshgrids
         * @param ddphi_xx array to store the second derivative wave functions in xx direction
         * @param ddphi_xy array to store the second derivative wave functions in xy direction
         * @param ddphi_xz array to store the second derivative wave functions in xz direction
         * @param ddphi_yy array to store the second derivative wave functions in yy direction
         * @param ddphi_yz array to store the second derivative wave functions in yz direction
         * @param ddphi_zz array to store the second derivative wave functions in zz direction
         */
        template <typename T>
        void set_ddphi(
            const std::vector<Vec3d>& coords, const int stride,
            T* ddphi_xx, T* ddphi_xy, T* ddphi_xz,
            T* ddphi_yy, T* ddphi_yz, T* ddphi_zz) const;

    private:
        // the atom object
        const Atom* atom_;

        // the global index of the atom
        int iat_;

        // the global index of the atom among the same type of atoms
        int ia_;

        // the index of big grid which contains this atom
        Vec3i biggrid_idx_;

        // the index of the unitcell which contains this atom
        Vec3i unitcell_idx_;

        // the relative Cartesian coordinates of this atom
        // with respect to the big grid that contains it
        Vec3d tau_in_biggrid_;

        // the numerical orbitals of this atom
        // In fact, I think the Numerical_Orbital class
        // should be a member of the Atom class, not the GintAtom class
        const Numerical_Orbital* orb_;

};

} // namespace ModuleGint