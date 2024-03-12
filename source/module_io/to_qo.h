#ifndef TOQO_H
#define TOQO_H

#include <iostream>
#include <string>
#include "module_cell/unitcell.h"
#include "module_basis/module_nao/two_center_integrator.h"
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
        3.1. calculate S(R) for all R: calculate_ovlpk(), calculate_ovlpR()
        3.2. fold S(R) to S(k): fold_ovlpR()
*/
class toQO
{
    public:
        // constructor is identical with what INPUT needs to define a toQO task
        // how user defines a QO task in INPUT, how it is constructed here, this
        // is to ensure the self-containedness of the class, as much as possible
        // to avoid introducing variables not from constructor or other under
        // control functions.
        toQO(const std::string& qo_basis,                   //< basis of QO, hydrogen or pswfc
             const std::vector<std::string>& strategies,    //< strategies for each atom type, more details see manual
             const double& qo_thr,                          //< threshold for QO
             const std::vector<double>& screening_coeffs);  //< screening coefficients for pseudowavefunction or Slater screening
        ~toQO();

        // initialize function is to import program-related information, it is,
        // more dynamic than constructor because this is actually an interface
        // to states of rest of program, unlike constructor, actually defines
        // the state of the class.
        void initialize(const std::string& out_dir,                                 //< directory of output files
                        const std::string& pseudo_dir,                              //< directory of pseudopotentials
                        const std::string& orbital_dir,                             //< directory of numerical atomic orbitals
                        const UnitCell* p_ucell,                                    //< interface to the unitcell
                        const std::vector<ModuleBase::Vector3<double>>& kvecs_d,    //< kpoints
                        std::ofstream& ofs_running,                                 //< output stream for running information
                        const int& rank,                                            //< rank of present processor
                        const int& nranks);                                         //< total number of processors
        // import structure, including supercell and kvectors are read in this function
        void read_structures(const UnitCell* p_ucell,                                   //< interface to the unitcell
                             const std::vector<ModuleBase::Vector3<double>>& kvecs_d,   //< kpoints
                             const int& iproc,                                          //< rank of present processor
                             const int& nprocs);                                        //< total number of processors
        
        // Two-center integral
        // QO is just one kind of representation, here it is representation from numerical 
        // atomic orbitals spanned space to atomic orbitals spanned space. Therefore two-center 
        // integral is the core of the whole transformation, it calculates the overlap between 
        // atomic orbitals and numerical atomic orbitals
        // for |>: to build RadialCollection filled with AtomicRadials
        void build_nao(const int ntype,                             //< number of atom types
                       const std::string orbital_dir,               //< directory of numerical atomic orbitals
                       const std::string* const orbital_fn,         //< filenames of numerical atomic orbitals
                       const int rank);                             //< rank of present processor
        // for <|: to build RadialCollection filled with PswfcRadials or HydrogenRadials
        void build_ao(const int ntype,                                                      //< number of atom types
                      const std::string pseudo_dir,                                         //< directory of pseudopotentials
                      const std::string* const pspot_fn = nullptr,                          //< filenames of pseudopotentials
                      const std::vector<double> screening_coeffs = std::vector<double>(),   //< screening coefficients of pseudopotentials
                      const double qo_thr = 1e-10,                                          //< threshold for QO
                      const std::ofstream& ofs = std::ofstream(),                           //< output stream for running information
                      const int rank = 0);                                                  //< rank of present processor
        // for <|: to build RadialCollection filled with HydrogenRadials
        void build_hydrogen(const int ntype,                    //< number of atom types
                            const double* const charges,        //< charges of atoms
                            const bool slater_screening,        //< whether use slater screening
                            const int* const nmax,              //< maximum principle quantum number of atoms
                            const double qo_thr,                //< threshold for QO
                            const int rank);                    //< rank of present processor
        // for <|: to build RadialCollection filled with PswfcRadials
        void build_pswfc(const int ntype,                           //< number of atom types
                         const std::string pseudo_dir,              //< directory of pseudopotentials
                         const std::string* const pspot_fn,         //< filenames of pseudopotentials
                         const double* const screening_coeffs,      //< screening coefficients of pseudopotentials, appears like a factor (exp[-s*r]) scaling the pswfc
                         const double qo_thr,                       //< threshold for QO
                         const int rank);                           //< rank of present processor
        // EXTERNAL EXPOSED FUNCTION, calculate all things in one shot
        void calculate();
        // calculate <A(i, R)|phi(j, R')> = Sij(R)
        void calculate_ovlpR(const int iR); //< iR index of supercell vector
        // calculate <A(i, k)|phi(j, k)> = Sij(k)
        void calculate_ovlpk(int ik);       //< ik index of kpoint
        // for overlap I/O
        // S to file
        template <typename T>
        void write_ovlp(const std::string& dir,             //< directory of output files
                        const std::vector<T>& matrix,       //< matrix to write
                        const int& nrows,                   //< number of rows
                        const int& ncols,                   //< number of columns
                        const bool& is_R = false,           //< whether it is in real space
                        const int& imat = 0);               //< index of matrix
        // S from file
        void read_ovlp(const std::string& dir,              //< directory of output files
                       const int& nrows,                    //< number of rows
                       const int& ncols,                    //< number of columns
                       const bool& is_R = false,            //< whether it is in real space
                       const int& imat = 0);                //< index of matrix
        /// @brief build bidirectional map indexing for one single RadialCollection object, which is an axis of two-center-integral table.
        /// @details from (it,ia,l,zeta,m) to index and vice versa
        void radialcollection_indexing(const RadialCollection&,                             //< [in] instance of RadialCollection
                                       const std::vector<int>&,                             //< [in] number of atoms for each type
                                       std::map<std::tuple<int,int,int,int,int>,int>&,      //< [out] mapping from (it,ia,l,zeta,m) to index
                                       std::map<int,std::tuple<int,int,int,int,int>>&);     //< [out] mapping from index to (it,ia,l,zeta,m)
        /// @brief calculate vectors connecting all atom pairs that needed to calculate their overlap
        ModuleBase::Vector3<double> cal_two_center_vector(ModuleBase::Vector3<double> rij,      //< vector connecting atom i and atom j
                                                          ModuleBase::Vector3<int> R);          //< supercell vector
        /// @brief for qo_basis pswfc only, return if include present orbitals in two-center integral
        bool orbital_filter(const int,              //< itype
                            const std::string);     //< notation of subshell, like all, s, p, d and any combination of them
        void deallocate_ovlp(const bool& is_R = false);     //< deallocate memory for ovlp_ao_nao_R_ or ovlp_ao_nao_k_
        void allocate_ovlp(const bool& is_R = false);       //< allocate memory for ovlp_ao_nao_R_ or ovlp_ao_nao_k_
        void zero_out_ovlps(const bool& is_R);              //< zero out ovlp_ao_nao_R_ or ovlp_ao_nao_k_
        void append_ovlpR_eiRk(int ik, int iR);             //< append S(R) to S(k), memory saving

        // MPI related
        void bcast_stdvector_ofvector3int(std::vector<ModuleBase::Vector3<int>>& vec);
        void bcast_stdvector_ofvector3double(std::vector<ModuleBase::Vector3<double>>& vec);

        // Neighboring list
        /// @brief get all possible (n1n2n3) defining supercell and scatter if MPI enabled
        void scan_supercell(const int& iproc,       //< rank of present processor
                            const int& nprocs);     //< total number of processors
        /// @brief this is a basic functional for scanning (ijR) pair for one certain i, return Rs
        /// @attention an algorithm loop over (i,)j,R, and return Rs
        /// @return a vector collects (n1, n2, n3) for present atom
        std::vector<ModuleBase::Vector3<int>> scan_supercell_for_atom(int it,               //< type of atom i
                                                                      int ia,               //< index of atom i
                                                                      int start_it = 0,     //< starting scan index of atom type
                                                                      int start_ia = 0);    //< starting scan index of atom
        /// @brief core algorithm to scan supercells, find the maximal supercell according to present cutoff radius
        /// @return a vector of (n1n2n3) defining supercell
        std::vector<int> rcut_to_supercell_index(double rcut,                       //< sum of cutoff radius of orbitals of atom i and atom j
                                                 ModuleBase::Vector3<double> a,     //< cell vector a (in Bohr)
                                                 ModuleBase::Vector3<double> b,     //< cell vector b (in Bohr)
                                                 ModuleBase::Vector3<double> c);    //< cell vector c (in Bohr)
        /// @brief get vector squared norm in supercell
        /// @return (rij + n1R1 + n2R2 + n3R3)^2
        double norm2_rij_supercell(ModuleBase::Vector3<double> rij,     //< vector connecting atom i and atom j in unitcell
                                   int n1,                              //< supercell index 1
                                   int n2,                              //< supercell index 2
                                   int n3);                             //< supercell index 3
        /// @brief eliminate duplicate vectors in a vector of vector3
        template <typename T>
        void eliminate_duplicate_vector3(std::vector<ModuleBase::Vector3<T>>& vector3s);
        /// @brief write supercells information to file
        void write_supercells();
        
        // getters
        int ntype() const { return ntype_; }
        int nks() const { return nks_; }
        std::string qo_basis() const { return qo_basis_; }
        std::vector<std::string> strategies() const { return strategies_; }
        std::string strategy(const int itype) const { return strategies_[itype]; }
        UnitCell* p_ucell() const { return const_cast<UnitCell*>(p_ucell_); }
        RadialCollection* p_nao() const { return nao_.get(); }
        RadialCollection* p_ao() const { return ao_.get(); }
        int nR() const { return nR_; }
        int nchi() const { return nchi_; }
        int nphi() const { return nphi_; }
        std::vector<ModuleBase::Vector3<int>> supercells() const { return supercells_; }
        std::vector<double> ovlpR() const { return ovlpR_; }
        double ovlpR(const int i, const int j) const { return ovlpR_[i*nchi_+j]; }
        std::vector<std::complex<double>> ovlpk() const { return ovlpk_; }
        std::complex<double> ovlpk(const int i, const int j) const { return ovlpk_[i*nchi_+j]; }
        std::vector<std::string> symbols() const { return symbols_; }
        std::vector<double> charges() const { return charges_; }
        atom_in atom_database() const { return atom_database_; }
        std::vector<ModuleBase::Vector3<double>> kvecs_d() const { return kvecs_d_; }

    private:
        // Variables defining QO task
        std::string qo_basis_ = "hydrogen";
        std::vector<std::string> strategies_;
        double qo_thr_ = 1e-10;
        std::vector<double> screening_coeffs_;
        // Variables defining I/O
        std::string out_dir_;                   //< directory of output files
        std::string pseudo_dir_;                //< directory of pseudopotentials
        std::string orbital_dir_;               //< directory of numerical atomic orbitals
        // Variables defining parallelism
        int iproc_ = 0;
        int nprocs_ = 1;
        // variables defining structure
        const UnitCell* p_ucell_ = nullptr;                        //< interface to the unitcell, its lifespan is not managed here
        std::vector<int> iRs_;                                     //< indices of supercell vectors (local)
        std::vector<ModuleBase::Vector3<int>> supercells_;         //< supercell vectors (global)
        std::vector<int> iks_;                                     //< indices of kpoints (local)
        std::vector<ModuleBase::Vector3<double>> kvecs_d_;         //< kpoints (global)
        // Two center integral
        std::unique_ptr<RadialCollection> nao_;                    //< numerical atomic orbitals
        std::unique_ptr<RadialCollection> ao_;                     //< atomic orbitals
        std::unique_ptr<TwoCenterIntegrator> overlap_calculator_;  //< two center integrator
        std::vector<double> ovlpR_;                                //< overlap between atomic orbitals and numerical atomic orbitals, in real space
        std::vector<std::complex<double>> ovlpk_;                  //< overlap between atomic orbitals and numerical atomic orbitals, in k space
        //< indices
        std::map<std::tuple<int,int,int,int,int>,int> index_ao_;   //< mapping from (it,ia,l,zeta,m) to index
        std::map<int,std::tuple<int,int,int,int,int>> rindex_ao_;  //< mapping from index to (it,ia,l,zeta,m)
        std::map<std::tuple<int,int,int,int,int>,int> index_nao_;  //< mapping from (it,ia,l,zeta,m) to index
        std::map<int,std::tuple<int,int,int,int,int>> rindex_nao_; //< mapping from index to (it,ia,l,zeta,m)
        // Variables defining dimensions or resource allocation
        int nks_ = 0;                                              //< number of kpoints for present processor, for S(k)
        int nks_tot_ = 0;                                          //< total number of kpoints
        int nR_ = 0;                                               //< number of supercell vectors on present processor, for S(R)
        int nR_tot_ = 0;                                           //< total number of supercell vectors
        int nchi_ = 0;                                             //< number of atomic orbitals, chi in \mathbf{S}^{\chi\phi}(\mathbf{k})
        int nphi_ = 0;                                             //< number of numerical atomic orbitals, phi in \mathbf{S}^{\chi\phi}(\mathbf{k})
        // Variables defining atoms
        atom_in atom_database_;                 //< atomic information database
        int ntype_ = 0;                         //< number of atom types
        std::vector<int> na_;                   //< number of atoms for each type
        std::vector<std::string> symbols_;      //< symbols of atoms
        std::vector<double> charges_;           //< charges of atoms
        std::vector<int> nmax_;                 //< maximum principle quantum number of atoms
};
#endif // TOQO_H