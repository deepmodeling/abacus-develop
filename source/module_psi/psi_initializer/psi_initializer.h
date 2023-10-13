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

template<typename T, typename Device>
class psi_initializer
{
    private:
        using Real = GetTypeReal<T>::type;
    public:
        psi_initializer(UnitCell* p_ucell_in,
                        ModulePW::PW_Basis_K* pw_wfc_in,
                        Parallel_Kpoints* p_parakpts_in,
                        Representation<T, Device>* p_rep_in): 
                        p_ucell(p_ucell_in), 
                        pw_wfc(pw_wfc_in),
                        p_parakpts(p_parakpts_in),
                        p_rep(p_rep_in) { };
        virtual ~psi_initializer();
        /*
            mainly used function (from external)
        */
        /// @brief initialize the given psi
        /// @param psi_out wavefunction to initialize
        /// @note initialize not fully initialize wavefunction, a very first diagonalization would be very beneficial
        virtual void initialize(const psi::Psi<T, Device>& psi_out) = 0;
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
        //#ifdef __MPI
        /// @brief (about planewaves distribution) from stick mapping to pool
        /// @param stick 
        /// @param ir 
        /// @param out 
        void stick_to_pool(Real* stick, const int& ir, Real* out) const;
        //#endif
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
        /// @brief setter of p_rep
        /// @param p_rep_in Representation pointer
        void set_interface_rep(Representation<T, Device>* p_rep_in) { this->p_rep = p_rep_in; }
        /// @brief getter of p_rep
        /// @return this->p_rep
        Representation<T, Device>* get_interface_rep() const { return this->p_rep; }
    protected:
        UnitCell* p_ucell;
        ModulePW::PW_Basis_K* pw_wfc;
        Parallel_Kpoints* p_parakpts;
        Representation<T, Device>* p_rep;
        /// @brief method of initializing psi
        std::string method = "random";
        /// @brief random_mix, the ratio of random number and original data
        double random_mix = 0.0;
    private:
        int nbands_complem = 0;
        int* ixy2is = nullptr;
        int random_seed = 0;
};