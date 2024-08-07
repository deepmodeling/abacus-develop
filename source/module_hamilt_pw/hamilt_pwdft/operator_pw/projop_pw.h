/**
 * A temporary demo for illustrating how to employ the real-space projection to implement
 * local orbital projector-based operators and derived algorithms.
 * 
 * This file will be removed after the implementation of DFT+U in PW basis.
 * 
 */
#include <string>
#include <vector>
#include "operator_pw.h"
#include "module_base/macros.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis_k.h"

namespace hamilt {

    #ifndef DFTUTEMPLATE_H
    #define DFTUTEMPLATE_H
    template<class T>
    class DFTU: public T
    {}; // this is a dummy class, it will be specialized later
    #endif
    template<typename T, typename Device>
    class DFTU<OperatorPW<T, Device>> : public OperatorPW<T, Device>
    {
        private:
            using Real = typename GetTypeReal<T>::type;
        public:
            // various constructors, support different types of projectors

            /**
             * @brief Construct a new DFTU object with interface similar with other
             * operator implementations
             * 
             * @param isk i am not clear what this is...
             * @param ucell_in pointer to unit cell, read-only accessible, lifespan not managed by this class
             * @param pw_basis pointer to PW_Basis_K, read-only accessible, lifespan not managed by this class
             */
            DFTU(const std::vector<int>& isk,
                 const UnitCell* ucell_in,      //< the dependency on UnitCell
                 const ModulePW::PW_Basis_K* pw_basis);

            /**
             * @brief Construct a new DFTU object with most flexible interface
             * 
             * @param isk i am not clear what this is...
             * @param l_hubbard angular momentum of the hubbard orbitals
             * @param u hubbard U values for each projector of each atom type
             * @param rgrid uniformed radial grid for the projectors
             * @param projs projectors for each atom type
             * @param natom number of atoms
             * @param tau atomic positions
             * @param omega cell volume
             * @param tpiba 2*pi/lat0
             * @param q q vector, {|k+G>}, may change its datatype for mult-k cases
             * @param dq for interpolation table, the interval of q
             * @param nq for interpolation table, the number of q-points
             */
            DFTU(const std::vector<int> isk,
                 const std::vector<int>& l_hubbard,
                 const std::vector<double>& u,
                 const std::vector<double>& rgrid,
                 const std::vector<std::vector<double>>& projs,
                 const std::vector<int>& natom,                         //< UnitCell::nat
                 const std::vector<ModuleBase::Vector3<double>*>& tau,  //< UnitCell::...
                 const double& omega,                                   //< UnitCell::omega
                 const double& tpiba,                                   //< UnitCell::tpiba
                 const std::vector<ModuleBase::Vector3<double>>& q,     //< PW_Basis::getgpluskcar
                 const double& dq = 0.01,                               //< GlobalV::DQ
                 const double& nq = 10000);                             //< GlobalV::NQX

            /**
             * @brief the CENTRAL function that can act on the wavefunction psi
             * 
             * @param nbands number of bands
             * @param nbasis number of basis
             * @param npol 
             * @param tmpsi_in input wavefunction
             * @param tmhpsi output/rotated wavefunction
             * @param ngk what??? number of |G> at present k point???
             */
            virtual void act(const int nbands,
                             const int nbasis,
                             const int npol,
                             const T* tmpsi_in,
                             T* tmhpsi,
                             const int ngk = 0) const override;

            /**
             * @brief read abacus numerical atomic orbital
             * 
             * @param forb file name
             * @param lmax the maximal angular momentum of orbitals defined in file
             * @param nzeta number of zeta for each angular momentum
             * @param nr number of radial grid points
             * @param dr radial grid spacing
             * @param radials radial orbitals
             */
            static void read_abacus_orb(const std::string& forb, 
                                        int& lmax, 
                                        std::vector<int>& nzeta,
                                        int& nr, 
                                        double& dr,
                                        std::vector<std::vector<double>>& radials);
        private:
            std::vector<std::complex<double>> proj_q_tab_; // stores the <beta|G+k>
    };
}