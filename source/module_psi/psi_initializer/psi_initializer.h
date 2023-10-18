// method support
#include "module_basis/module_representation/representation.h"
#include "module_hsolver/diago_iter_assist.h"
// data support
#include "module_base/global_variable.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
// data structure support
#include "module_basis/module_pw/pw_basis_k.h"
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
// template programming support
#include "module_base/macros.h"

/*
Psi (planewave based wavefunction) initializer
Auther: Kirk0830
Institute: AI for Science Institute, BEIJING

This class is used to allocate memory and give initial guess for psi (not kspw_psi the FPTYPE, Device template one)
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

How does this class work?
The central function is initialize(), 
For init_wfc = random, initialize wavefunction directly
*/

template<typename T, typename Device>
class psi_initializer
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        #ifdef __MPI
        psi_initializer(UnitCell* p_ucell_in,
                        ModulePW::PW_Basis_K* pw_wfc_in,
                        Parallel_Kpoints* p_parakpts_in,
                        Representation<T, Device>* p_rep_in);
        #else
        psi_initializer(UnitCell* p_ucell_in,
                        ModulePW::PW_Basis_K* pw_wfc_in,
                        Representation<T, Device>* p_rep_in);
        #endif
        ~psi_initializer();
        /*
            mainly used function (from external)
        */
        /// @brief initialize the given psi
        /// @param psi_out wavefunction to initialize
        /// @note initialize not fully initialize wavefunction, a very first diagonalization would be very beneficial
        void initialize();
        /// @brief allocate memory for both psi::Psi<std::complex<double>, Device> ESolver_FP::psi and Representation<T, Device>::psig
        /// @return pointer/memory address of allocated ESolver_FP::psi
        psi::Psi<std::complex<double>>* allocate();
        /*
            used function (internal)
        */
        // random to complement bands not initialized by pswfc or nao, therefore it is a basic function, or psi_initializer_random will be inherented by all other methods.
        /// @brief kernel to generate and assign random number for psi
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        void random_t(T* psi, const int iw_start, const int iw_end, const int ik);
        #ifdef __MPI
        /// @brief (about planewaves distribution) from stick mapping to pool
        /// @param stick 
        /// @param ir 
        /// @param out 
        void stick_to_pool(Real* stick, const int& ir, Real* out) const;
        #endif
        /*
            setters and getters
        */
        /// @brief get complement number of bands
        /// @return nbands_complem
        int get_nbands_complem() const { return this->nbands_complem; }
        /// @brief set number of complementary bands
        /// @param nbands_in nbands_complem
        void set_nbands_complem(int nbands_in) { this->nbands_complem = nbands_in; }
        /// @brief get method of initializing psi
        /// @return the method
        std::string get_method() const { return this->method; }
        /// @brief set method manually
        /// @param method_in initialization method
        void set_method(std::string method_in) { this->method = method_in; }
        /// @brief setter of random_mix
        /// @param random_mix_in new value of random_mix
        void set_random_mix(const double random_mix_in) { this->random_mix = random_mix_in; }
        /// @brief getter of random_mix
        /// @return this->random_mix
        double get_random_mix() const { return this->random_mix; }
        /// @brief getter of ixy2is, the mapping from fftixy to stick index
        /// @return this->ixy2is
        int* get_ixy2is() const { return this->ixy2is; }
        /// @brief setter of ixy2is, the mapping from fftixy to stick index
        void set_ixy2is(int* ixy2is_in) { this->ixy2is = ixy2is_in; }
        /// @brief setter of random_seed
        /// @param random_seed_in new value of random_seed
        void set_random_seed(const int random_seed_in) { this->random_seed = random_seed_in; }
        /// @brief getter of random_seed
        /// @return this->random_seed
        int get_random_seed() const { return this->random_seed; }
        /// @brief setter of p_ucell
        /// @param p_ucell_in UnitCell pointer
        void set_interface_ucell(UnitCell* p_ucell_in) { this->p_ucell = p_ucell_in; }
        /// @brief getter of p_ucell
        /// @return this->p_ucell
        UnitCell* get_interface_ucell() const { return this->p_ucell; }
        /// @brief setter of pw_wfc
        /// @param pw_wfc_in PW_Basis_K pointer
        void set_interface_pw_wfc(ModulePW::PW_Basis_K* pw_wfc_in) { this->pw_wfc = pw_wfc_in; }
        /// @brief getter of pw_wfc
        /// @return this->pw_wfc
        ModulePW::PW_Basis_K* get_interface_pw_wfc() const { return this->pw_wfc; }
        #ifdef __MPI
        /// @brief setter of p_parakpts
        /// @param p_parakpts_in Parallel_Kpoints pointer
        void set_interface_parakpts(Parallel_Kpoints* p_parakpts_in) { this->p_parakpts = p_parakpts_in; }
        /// @brief getter of p_parakpts
        /// @return this->p_parakpts
        Parallel_Kpoints* get_interface_parakpts() const { return this->p_parakpts; }
        #endif
        /// @brief setter of p_rep
        /// @param p_rep_in Representation pointer
        void set_interface_rep(Representation<T, Device>* p_rep_in) { this->p_rep = p_rep_in; }
        /// @brief getter of p_rep
        /// @return this->p_rep
        Representation<T, Device>* get_interface_rep() const { return this->p_rep; }
    protected:
        UnitCell* p_ucell = nullptr;
        ModulePW::PW_Basis_K* pw_wfc = nullptr;
        #ifdef __MPI
        Parallel_Kpoints* p_parakpts = nullptr;
        #endif
        Representation<T, Device>* p_rep = nullptr;
        /// @brief method of initializing psi
        std::string method = "random";
        /// @brief random_mix, the ratio of random number and original data
        Real random_mix = 0.05;
    private:
        int mem_saver = 0; // will deprecated this variable soon
        int nbands_complem = 0;
        int* ixy2is = nullptr;
        int random_seed = 0;
};