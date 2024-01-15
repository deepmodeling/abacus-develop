#ifndef TOQO_H
#define TOQO_H

#include <iostream>
#include <string>
#include "module_cell/unitcell.h"
#include "module_basis/module_nao/two_center_integrator.h"
#include "module_basis/module_pw/pw_basis_k.h"
#include "module_base/atom_in.h"
#include "module_base/vector3.h"
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

    Workflow:
    1. create an instance of toQO: toQO()
    2. initialize the class: initialize()
        2.1. import information from unitcell: unwrap_unitcell()
        2.2. build orbitals, nao and ao: build_nao() and build_ao()
        2.3. scan neighboring list: scan_supercell()
        2.4. allocate memory for ovlp_ao_nao_R_ and ovlp_ao_nao_k_: clean_up(), deallocate_ovlp() and allocate_ovlp()
    3. calculate overlap S(k): calculate()
        3.1. calculate S(R) for all R: calculate_ovlp_k(), calculate_ovlp_R()
        3.2. fold S(R) to S(k): fold_ovlp_R()
*/
class toQO
{
    private:
        using CplxMatrix = std::vector<std::vector<std::complex<double>>>;
        using RealMatrix = std::vector<std::vector<double>>;

    public:
        toQO(std::string qo_basis, std::string strategy = "minimal");
        ~toQO();

        /*
         *  Initialization
        */
        /// @brief initialize the QO class
        /// @param p_ucell interface (raw pointer) to the unitcell
        /// @param nkpts number of kpoints
        void initialize(UnitCell* p_ucell,
                        const std::vector<ModuleBase::Vector3<double>>& kvecs_d);
        /// @brief to get rid of direct use of UnitCell
        /// @param p_ucell 
        void unwrap_unitcell(UnitCell* p_ucell);
        /*
         *  Two center integrator interfaces
        */
        /// @brief build RadialCollection for numerical atomic orbitals
        /// @param ntype number of atom types
        /// @param orbital_fn filenames of numerical atomic orbitals
        void build_nao(const int ntype, const std::string* const orbital_fn);
        /// @brief build RadialCollection for atomic orbitals
        /// @param ntype number of atom types
        /// @param charges charges of atoms
        /// @param nmax maximum principle quantum number of atoms
        void build_hydrogen(const int ntype, const double* const charges, const int* const nmax);
        /// @brief build RadialCollection for atomic orbitals
        /// @param ntype number of atom types
        /// @param pspot_fn filenames of pseudopotentials
        /// @param screening_coeffs screening coefficients of pseudopotentials, appears like a factor (exp[-s*r]) scaling the pswfc
        void build_pswfc(const int ntype, const std::string* const pspot_fn, const double* const screening_coeffs);
        /// @brief build RadialCollection for atomic orbitals
        /// @param ntype number of atom types
        /// @param pspot_fn filenames of pseudopotentials, if use qo_basis = hydrogen, omit this parameter
        void build_ao(const int ntype, const std::string* const pspot_fn = nullptr);
        /*
         *   Main functions
        */
        /// @brief calculate the overlap between atomic orbitals and numerical atomic orbitals, in real space, at R[iR]
        /// @param iR index of supercell vector
        /// @note to save memory, the workflow can be organized as, once one S(R) is calculated, fold it to S(k), then clean up S(R)...
        void calculate_ovlp_R(const int iR);
        /// @brief calculate the overlap between atomic orbitals and numerical atomic orbitals, in k space
        /// @param kvec_c vector3 specifying a kpoint
        void calculate_ovlp_k(int ik);
        /// @brief calculate the overlap between atomic orbitals and numerical atomic orbitals, in k space and write to file
        void calculate();
        /// @brief write two dimensional matrix to file
        /// @tparam T type of matrix
        /// @param matrix matrix to write
        /// @param ik index of kpoint
        template <typename T>
        void write_ovlp(const std::vector<std::vector<T>>& matrix, const int& ik);
        /// @brief write supercells information to file
        void write_supercells();
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
        /// @brief calculate vectors connecting all atom pairs that needed to calculate their overlap
        ModuleBase::Vector3<double> cal_two_center_vector(ModuleBase::Vector3<double> rij,
                                                          ModuleBase::Vector3<int> R);
        /// @brief this is a basic functional for scanning (ijR) pair for one certain i, return Rs
        /// @attention an algorithm loop over (i,)j,R
        /// @param it type of atom i
        /// @param ia index of atom i
        /// @param start_it starting scan index of atom type
        /// @param start_ia starting scan index of atom
        /// @param rcut cutoff radius of numerical atomic orbital of atom i
        /// @return a vector collects (n1, n2, n3) for present atom
        std::vector<ModuleBase::Vector3<int>> scan_supercell_for_atom(int it, int ia, int start_it = 0, int start_ia = 0);
        /// @brief core algorithm to scan supercells, find the maximal supercell according to present cutoff radius
        /// @param rcut sum of cutoff radius of orbitals of atom i and atom j
        /// @param a cell vector a (in Bohr)
        /// @param b cell vector b (in Bohr)
        /// @param c cell vector c (in Bohr)
        /// @return a vector of (n1n2n3) defining supercell
        std::vector<int> rcut_to_supercell_index(double rcut, ModuleBase::Vector3<double> a, ModuleBase::Vector3<double> b, ModuleBase::Vector3<double> c);
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
        /// @brief eliminate duplicate vectors in a vector of vector3
        /// @tparam T type of vector3
        /// @param vector3s vector of vector3, both input and output
        template <typename T>
        void eliminate_duplicate_vector3(std::vector<ModuleBase::Vector3<T>>& vector3s);
        /*
            Data management
        */
        /// @brief clean up ovlp_ao_nao_R_ or ovlp_ao_nao_k_
        /// @param is_R true for ovlp_ao_nao_R_, false for ovlp_ao_nao_k_
        void deallocate_ovlp(bool is_R = false);
        /// @brief allocate memory for ovlp_ao_nao_R_ or ovlp_ao_nao_k_
        /// @param is_R true for ovlp_ao_nao_R_, false for ovlp_ao_nao_k_
        void allocate_ovlp(bool is_R = false);
        /// @brief clean up all data members
        void clean_up();
        /// @brief zero out ovlp_ao_nao_R_ or ovlp_ao_nao_k_
        /// @param is_R true for ovlp_ao_nao_R_, false for ovlp_ao_nao_k_
        void zero_out_ovlps(const bool is_R);

        /// @brief given a vector3 specifying a kpoint, fold ovlp_ao_nao_R_ (series of S(R), memory consuming)
        /// @param ik index of vector3 specifying a kpoint
        void fold_ovlp_R(int ik);
        /// @brief given a vector3 specifying a kpoint, append one single S(R), multiply by exp(-i*k*R) and add to ovlp_ao_nao_k_
        /// @param ik index of  vector3 specifying a kpoint
        /// @param iR index of supercell vector
        void append_ovlp_R_eiRk(int ik, int iR);

        // setters
        void set_qo_basis(const std::string qo_basis) { qo_basis_ = qo_basis; }
        void set_strategy(const std::string strategy) { strategy_ = strategy; }
        void set_save_mem(const bool save_mem) { save_mem_ = save_mem; }
        
        // getters
        int ntype() const { return ntype_; }
        int nkpts() const { return nkpts_; }
        std::string qo_basis() const { return qo_basis_; }
        std::string strategy() const { return strategy_; }
        UnitCell* p_ucell() const { return p_ucell_; }
        RadialCollection* p_nao() const { return nao_.get(); }
        RadialCollection* p_ao() const { return ao_.get(); }
        int nR() const { return nR_; }
        int nchi() const { return nchi_; }
        int nphi() const { return nphi_; }
        std::vector<ModuleBase::Vector3<int>> supercells() const { return supercells_; }
        std::vector<RealMatrix> ovlp_R() const { return ovlp_R_; }
        CplxMatrix ovlp_k() const { return ovlp_k_; }
        bool save_mem() const { return save_mem_; }
        std::vector<std::string> symbols() const { return symbols_; }
        std::vector<double> charges() const { return charges_; }
        atom_in atom_database() const { return atom_database_; }
        std::vector<ModuleBase::Vector3<double>> kvecs_d() const { return kvecs_d_; }

    private:
        //
        // interface to deprecated
        //
        /// @brief interface to the unitcell
        UnitCell* p_ucell_ = nullptr;

        //
        // high dimensional data
        //
        /// @brief supercell vectors
        std::vector<ModuleBase::Vector3<int>> supercells_;
        /// @brief overlaps between atomic orbitals and numerical atomic orbitals, in real space
        std::vector<RealMatrix> ovlp_R_;
        /// @brief overlap between atomic orbitals and numerical atomic orbitals, in k space
        CplxMatrix ovlp_k_;

        //
        // basic data member
        //
        /// @brief number of kpoints, for S(k)
        int nkpts_ = 0;
        /// @brief number of supercell vectors, for S(R)
        int nR_ = 0;
        /// @brief number of atomic orbitals, chi in \mathbf{S}^{\chi\phi}(\mathbf{k})
        int nchi_ = 0;
        /// @brief number of numerical atomic orbitals, phi in \mathbf{S}^{\chi\phi}(\mathbf{k})
        int nphi_ = 0;
        //
        // data unwrapped from unitcell
        //
        /// @brief number of atom types
        int ntype_ = 0;
        /// @brief number of atoms for each type
        std::vector<int> na_;
        //
        // qo_basis = strategy
        //
        std::vector<std::string> symbols_;
        std::vector<double> charges_;
        std::vector<int> nmax_;
        //
        //
        //
        std::vector<ModuleBase::Vector3<double>> kvecs_d_;
        //
        // attributes
        //
        /// @brief current atomic orbital basis for generating QO
        /// @details hydrogen: hydrogen-like orbitals; pswfc: pseudowavefunction
        std::string qo_basis_ = "hydrogen";
        /// @brief strategy for generating QO
        /// @details full: 1s, 2s, 2p, 3s, 3p, 3d, ...
        ///          minimal: 1s, 2p, 3d, 4f, ...
        ///          energy: according to Hund's rule
        std::string strategy_ = "minimal";

        //
        // memory control
        //
        /// @brief whether to save memory
        bool save_mem_ = true;

        // 
        // advanced data structures
        //
        // orbitals data carriers
        /// @brief numerical atomic orbitals
        std::unique_ptr<RadialCollection> nao_;
        /// @brief atomic orbitals
        std::unique_ptr<RadialCollection> ao_;
        /// @brief two center integrator
        std::unique_ptr<TwoCenterIntegrator> overlap_calculator_;
        /// @brief atom database
        atom_in atom_database_;
};
#endif // TOQO_H