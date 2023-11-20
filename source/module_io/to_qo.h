#ifndef TOQO_H
#define TOQO_H

#include <iostream>
#include <string>
#include "module_cell/unitcell.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_base/atom_in.h"
/*
    Quasiatomic Orbital (QO) transformation and analysis

    Technical details and procedures:
    1. first project a given set of atomic orbitals onto one certain set, nao or pw(not implemented)
       AO now supports hydrogen-like radial orbitals
    2. project AO onto the eigenstate, then substrate them from AO, get AO'
    3. Diagonalize AO' and use canoincal orthogonalization to find in total m states (m is arbitrary to user)
    4. merge space spanned by eigenstate and the one of m states
    5. project AO onto this new space basis

    Functionalities:
    1. support different projectors like:
      (1) hydrogen-like: full, minimal and double. Full may exceed the dimensional of space of nao, minimal
          will not in most cases
      (2) pseudowavefunction (pswfc): only avaiable for pseudpotentials with pswfc, unlike SG15.
    2. output overlap between the given projectors and nao, so make it easy to calculate QO and Hamiltonian in
      QO representation.
*/
class toQO
{
    private:
        using Matrix = std::vector<std::vector<std::complex<double>>>;

    public:
        toQO(std::string qo_basis, std::string strategy = "");
        ~toQO();
        // ----
        // main functions, implemented in to_qo.cpp
        // ----
        /// @brief initialize the QO class
        /// @param p_ucell interface (raw pointer) to the unitcell
        /// @param nkpts number of kpoints
        void initialize(UnitCell *p_ucell,
                        ModulePW::PW_Basis_K* p_pw_wfc,
                        int nkpts);
        /// @brief calculate the overlap between atomic orbitals and numerical atomic orbitals, in real space
        void cal_ovlp_ao_nao_R(const int iR);
        /// @brief calculate the overlap between atomic orbitals and numerical atomic orbitals, in k space
        void cal_ovlp_ao_nao_k(const int ik);
        /// @brief write the overlap to file
        void write_ovlp_ao_nao();

        // ----
        // tool functions, implemented in to_qo_tools.cpp
        // ----
        /// @brief calculate vectors connecting all atom pairs that needed to calculate their overlap
        void cal_two_center_vectors();

        /*
            Neighboring list searching algorithm (not implemented yet)

            The neighboring list here is for searching all possible (ijR) pairs that overlap between atom i and atom j in different
            cells distanced by R still has nonzero value. Therefore it is actually a problem:

            |rij + n1R1 + n2R2 + n3R3| >= (rcut,i + rcut,j),
            , where n1, n2, n3 are integers, R1, R2, R3 are lattice vectors, and rcut,i, rcut,j are cutoff radii of numerical orbitals 
            of atom i and atom j. rij is the distance between atom i and atom j.in the unitcell.
            Take square on both sides, we have
            rij^2 + 2 rij * (n1R1 + n2R2 + n3R3) + (n1R1 + n2R2 + n3R3)^2 >= (rcut,i + rcut,j)^2
            . n1, n2 and n3 are values of interest, rij and rcut,i and rcut,j are known for specific atom pair ij.

            Therefore neighboring list searching problem is a problem of root finding of a quadratic equation.
            The quadratic equation is
            (R1^2)*n1^2 + (R2^2)*n2^2 + (R3^2)*n3^2 + 2*(R1*R2)*n1*n2 + 2*(R1*R3)*n1*n3 + 2*(R2*R3)*n2*n3
            + 2*rij*(R1*n1 + R2*n2 + R3*n3) + rij^2 - (rcut,i + rcut,j)^2 = 0
            one can drop the constraint that n1, n2 and n3 are integers, then use ceiling to get the integer values.

            To solve the quadratic equation, one can rotate the coordinate system so that the function can become
            a sphere. Then it will be an approximately 2d scan (phi, theta) instead of a 3d scan (n1, n2, n3). The time complexity will be reduced from 
            O(Natom^2*Ncell^3) to O(Natom^2*Ncell^2).

            This algorithm is not implemented yet.
            A diagonalization of third order matrix of coefficients of ninj can transform the targeting function
            to a translated sphere, then translate it back to origin. Then it will be a 2d scan (phi, theta).
        */
        /// @brief this is a basic functional for scanning (ijR) pair for one certain i, return Rs
        /// @attention an algorithm loop over (i,)j,R
        /// @param it type of atom i
        /// @param ia index of atom i
        /// @param start_it starting scan index of atom type
        /// @param start_ia starting scan index of atom
        /// @param rcut cutoff radius of numerical atomic orbital of atom i
        /// @return a vector collects (n1, n2, n3) for present atom
        std::vector<ModuleBase::Vector3<int>> scan_supercell_for_atom(int it, int ia, int start_it = 0, int start_ia = 0);

        /// @brief get vector squared norm in supercell
        /// @param rij rij in unitcell
        /// @param n1 supercell index 1
        /// @param n2 supercell index 2
        /// @param n3 supercell index 3
        /// @return (rij + n1R1 + n2R2 + n3R3)^2
        double norm2_rij_supercell(ModuleBase::Vector3<double> rij, int n1, int n2, int n3);
        
        /// @brief get all possible (n1n2n3) defining supercell
        /// @return a vector of (n1n2n3)
        void scan_supercell();

        Matrix folding_matrixR(ModuleBase::Vector3<double> kvec_c);

    private:

        /// @brief interface to the unitcell
        UnitCell *p_ucell_ = nullptr;
        /// @brief interface to the pw wavefunction
        ModulePW::PW_Basis_K *p_pw_wfc_ = nullptr;
        /// @brief number of kpoints
        int nkpts_ = 0;

        std::vector<ModuleBase::Vector3<int>> supercells_;

        std::vector<Matrix> ovlp_ao_nao_R_;
        std::vector<Matrix> ovlp_ao_nao_k_;

        /// @brief current atomic orbital basis for generating QO
        /// @details hydrogen_minimal: 1s, 2p, 3d, ... 
        ///          hydrogen_double: 1s, 2s, 2p, 3p, 3d, ...
        ///          hydrogen: 1s, 2s, 2p, 3s, 3p, 3d, ... 
        std::string qo_basis_ = "hydrogen";
        std::string strategy_ = "minimal";

        std::unique_ptr<RadialCollection> nao_;
        std::unique_ptr<RadialCollection> ao_;

        atom_in atom_database_;
};
#endif // TOQO_H