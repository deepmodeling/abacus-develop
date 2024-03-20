#ifndef PSI_INITIALIZER_H
#define PSI_INITIALIZER_H
// data structure support
#include "module_psi/psi.h" // for psi data structure
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
// smart pointer for auto-memory management
#include <memory>
// numerical algorithm support
#include "module_base/spherical_bessel_transformer.h" // for spherical bessel transform
#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#endif

#include "module_base/macros.h"
#include <type_traits>
/*
Psi (planewave based wavefunction) initializer
Auther: Kirk0830
Institute: AI for Science Institute, BEIJING

This class is used to allocate memory and give initial guess for psi
therefore only double datatype is needed to be supported.
Following methods are available:
    1. random: use random number to initialize psi
               implemented in psi_initializer_random.h
    2. atomic: use pseudo-wavefunction in pseudopotential file to initialize psi
               implemented in psi_initializer_atomic.h
    3. atomic+random: mix 'atomic' with some random numbers to initialize psi
    4. nao: use numerical orbitals to initialize psi
            implemented in psi_initializer_nao.h
    5. nao+random: mix 'nao' with some random numbers to initialize psi
*/
template<typename T, typename Device>
class psi_initializer
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        // technical notes:
        // Polymorphism is used to implement different methods, and achieved by pointers and virtual functions
        psi_initializer() {};
        virtual ~psi_initializer() {};
        #ifdef __MPI // MPI additional implementation
        /// @brief initialize the psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,              //< structure factor
                                ModulePW::PW_Basis_K*,          //< planewave basis
                                UnitCell*,                      //< unit cell
                                Parallel_Kpoints* = nullptr,    //< parallel kpoints
                                const int& = 1,                 //< random seed
                                pseudopot_cell_vnl* = nullptr,  //< nonlocal pseudopotential
                                const int& = 0) = 0;            //< MPI rank
        
        Parallel_Kpoints* p_parakpts() const { return this->p_parakpts_; }
        void set_parakpts(Parallel_Kpoints* p_parakpts) { this->p_parakpts_ = p_parakpts; }
        /// @brief mapping from (ix, iy) to is
        void stick_to_pool(Real* stick,         //< stick
                           const int& ir,       //< ir
                           Real* out) const;    //< out
        #else
        /// @brief serial version of initialize function, link psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,                  //< structure factor
                                ModulePW::PW_Basis_K*,              //< planewave basis
                                UnitCell*,                          //< unit cell
                                const int& = 1,                     //< random seed
                                pseudopot_cell_vnl* = nullptr) = 0; //< nonlocal pseudopotential
        #endif
        /// @brief CENTRAL FUNCTION: allocate memory for psi
        psi::Psi<std::complex<double>>* allocate(bool only_psig = false);           //< if only allocate memory for psig

        void random_t(T* psi,                       //< [out] psi
                      const int iw_start,           //< iw_start, starting band index
                      const int iw_end,             //< iw_end, ending band index
                      const int ik);                //< ik, kpoint index
        virtual void random(T* psi,                 //< [out] psi
                            const int iw_start,     //< iw_start, starting band index
                            const int iw_end,       //< iw_end, ending band index
                            const int ik)           //< ik, kpoint index
        { ModuleBase::WARNING_QUIT("psi_initializer::random", "Polymorphism error"); }
        /// @brief CENTRAL FUNCTION: import external files, such as pseudopotential file or numerical orbitals file
        virtual void read_external_orbs(std::string* external_files,    //< orbital files
                                        const int& rank = 0)            //< MPI rank
        { ModuleBase::WARNING_QUIT("psi_initializer::set_orbital_files", "Polymorphism error"); }
        /// @brief CENTRAL FUNCTION: allocate interpolate table recording overlap integral between radial function and Spherical Bessel function
        virtual void allocate_table() { ModuleBase::WARNING_QUIT("psi_initializer::create_ovlp_table", "Polymorphism error"); }
        /// @brief CENTRAL FUNCTION: calculate the interpolate table
        virtual void tabulate() = 0;
        /// @brief CENTRAL FUNCTION: calculate projection of atomic radial function onto planewave basis BASED ON THE OVERLAP TABLE
        virtual void proj_ao_onkG(int ik) = 0;

        // getter and setter
        UnitCell* p_ucell() const { return this->p_ucell_; }
        pseudopot_cell_vnl* p_pspot_nl() const { return this->p_pspot_nl_; }
        Structure_Factor* p_sf() const { return this->sf_; }
        ModulePW::PW_Basis_K* pw_wfc() const { return this->pw_wfc_; }
        int random_seed() const { return this->random_seed_; }
        std::vector<int> ixy2is() const { return this->ixy2is_; }
        int mem_saver() const { return this->mem_saver_; }
        double random_mix() const { return this->random_mix_; }
        bool initialized() const { return this->initialized_; }
        std::string method() const { return this->method_; }
        int nbands_complem() const { return this->nbands_complem_; }
        std::weak_ptr<psi::Psi<T, Device>> share_psig() { return this->psig_; }

        void set_ucell(UnitCell* p_ucell_in) { this->p_ucell_ = p_ucell_in; }
        void set_pspot_nl(pseudopot_cell_vnl* p_pspot_nl_in) { this->p_pspot_nl_ = p_pspot_nl_in; }
        void set_sf(Structure_Factor* sf_in) { this->sf_ = sf_in; }
        void set_pw_wfc(ModulePW::PW_Basis_K* pw_wfc_in) { this->pw_wfc_ = pw_wfc_in; }
        void set_random_mix(const double random_mix_in) { this->random_mix_ = random_mix_in; }
        void set_ixy2is(const std::vector<int>& ixy2is_in) { this->ixy2is_ = ixy2is_in; }
        void set_random_seed(const int random_seed_in) { this->random_seed_ = random_seed_in; }
        void set_mem_saver(const int mem_saver_in) { this->mem_saver_ = mem_saver_in; }
        void set_initialized(bool initialized_in) { this->initialized_ = initialized_in; }
        void set_method(std::string method_in) { this->method_ = method_in; }
        void set_nbands_complem(int nbands_in) { this->nbands_complem_ = nbands_in; }

        // tool methods
        template <typename U>
        typename std::enable_if<std::is_same<U, float>::value, U>::type cast_to_T(const std::complex<double> in) {return static_cast<float>(in.real());}
        template <typename U>
        typename std::enable_if<std::is_same<U, double>::value, U>::type cast_to_T(const std::complex<double> in) {return static_cast<double>(in.real());}
        template <typename U>
        typename std::enable_if<std::is_same<U, std::complex<float>>::value, U>::type cast_to_T(const std::complex<double> in) {return std::complex<float>(static_cast<float>(in.real()), static_cast<float>(in.imag()));}
        template <typename U>
        typename std::enable_if<std::is_same<U, std::complex<double>>::value, U>::type cast_to_T(const std::complex<double> in) {return std::complex<double>(in.real(), in.imag());}

    protected:
        Structure_Factor* sf_ = nullptr;
        ModulePW::PW_Basis_K* pw_wfc_ = nullptr;
        UnitCell* p_ucell_ = nullptr;
        #ifdef __MPI
        Parallel_Kpoints* p_parakpts_ = nullptr;
        #endif
        pseudopot_cell_vnl* p_pspot_nl_ = nullptr;
        // shared by atomic, nao, atomic+random, nao+random
        ModuleBase::SphericalBesselTransformer sbt; // useful for atomic-like methods
        // shared by random, atomic+random, nao+random
        int random_seed_ = 1;
        // in old version it is of datatype int*, use std::vector<int> to avoid memory leak
        std::vector<int> ixy2is_;

        std::shared_ptr<psi::Psi<T, Device>> psig_;
    private:
        int mem_saver_ = 0;
        std::string method_ = "none";
        int nbands_complem_ = 0;
        bool initialized_ = false;
        double random_mix_ = 0;
};
#endif