#pragma once

#include "module_cell/atom_spec.h"
#include "module_basis/module_ao/ORB_atomic.h"
#include "gint_type.h"

namespace Gint
{

class GintAtom
{
    public:
        // constructor
        GintAtom(
            const Atom* atom,
            int iat,
            Vec3i biggrid_idx,
            Vec3i unitcell_idx,
            Vec3d tau_in_biggrid,
            const Numerical_Orbital* orb)
            : atom_(atom), iat_(iat), biggrid_idx_(biggrid_idx),
              unitcell_idx_(unitcell_idx), tau_in_biggrid_(tau_in_biggrid),
              orb_(orb){};

        // getter functions
        const Atom* get_atom() const { return atom_; };
        const int get_iat() const { return iat_; };
        const Vec3i& get_biggrid_idx() const { return biggrid_idx_; };
        const Vec3i& get_unitcell_idx() const { return unitcell_idx_; };
        const Vec3i& get_r() const { return unitcell_idx_; };
        const Vec3d& get_tau_in_biggrid() const { return tau_in_biggrid_; };
        const Numerical_Orbital* get_orb() const { return orb_; };

        const int get_nw() const { return atom_->nw; };
        const double get_rcut() const { return orb_->getRcut(); };
        
        /**
         * @brief Get the wave function values of the atom at a meshgrid.
         * 
         * psi[(n-1)*stride] ~ psi[(n-1)*stride + nw] store the wave function values of the first atom at the nth meshgrid
         * 
         * @param coords the cartesian coordinates of the meshgrids of a biggrid relative to the atom
         * @param stride the stride of the psi array between two adjacent meshgrids
         * @param psi array to store the wave function values
         */
        template <typename T>
        void get_psi(const std::vector<Vec3d> coords, const int stride, T* psi);

        /**
         * @brief Get the wave function values and its derivative
         * 
         * The reason for combining the functions to solve the wave function values 
         * and wave function derivatives into one function is to improve efficiency.
         * psi[(n-1)*stride] ~ psi[(n-1)*stride + nw] store the wave function values of the first atom at the nth meshgrid
         * 
         * @param coords the cartesian coordinates of the meshgrids of a biggrid relative to the atom
         * @param stride the stride of the psi array between two adjacent meshgrids
         * @param psi array to store the wave function values
         * @param dpsi_x array to store the derivative wave functions in x direction
         * @param dpsi_y array to store the derivative wave functions in y direction
         * @param dpsi_z array to store the derivative wave functions in z direction
         */
        template <typename T>
        void get_psi_dpsir(
            const std::vector<Vec3d> coords, const int stride,
            T* psi, T* dpsi_x, T* dpsi_y, T* dpsi_z);

        /**
         * @brief Get the wave function values and its second derivative
         * 
         * ddpsi[(n-1)*stride] ~ ddpsi[(n-1)*stride + nw] store the second derivative of 
         * wave function values of the atom at the first meshgrid
         *  
         * @param coords the cartesian coordinates of the meshgrids of a biggrid relative to the atom
         * @param stride the stride of the psi array between two adjacent meshgrids
         * @param ddpsi_xx array to store the second derivative wave functions in xx direction
         * @param ddpsi_xy array to store the second derivative wave functions in xy direction
         * @param ddpsi_xz array to store the second derivative wave functions in xz direction
         * @param ddpsi_yy array to store the second derivative wave functions in yy direction
         * @param ddpsi_yz array to store the second derivative wave functions in yz direction
         * @param ddpsi_zz array to store the second derivative wave functions in zz direction
         */
        template <typename T>
        void get_ddpsir(
            const std::vector<Vec3d> coords, const int stride,
            T* ddpsi_xx, T* ddpsi_xy, T* ddpsi_xz,
            T* ddpsi_yy, T* ddpsi_yz, T* ddpsi_zz);

    private:
        // the atom object
        const Atom* atom_;

        // the global index of the atom
        int iat_;

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

}

} // namespace Gint