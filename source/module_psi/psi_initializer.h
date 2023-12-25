#ifndef PSI_INITIALIZER_H
#define PSI_INITIALIZER_H
// data structure support
#include "module_psi/psi.h" // for psi data structure
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
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
*/
template<typename T, typename Device>
class psi_initializer
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        /// @brief default constructor of psi initializer
        psi_initializer() { };
        #ifdef __MPI
        /// @brief parameterized constructor of psi initializer (with MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param p_parakpts_in interface, link with Parallel_Kpoints GlobalC::Pkpoints
        /// @param random_seed_in random seed
        psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in = 1);
        #else
        /// @brief parameterized constructor of psi initializer (without MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param random_seed_in random seed
        psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in = 1);
        #endif
        /// @brief destructor
        virtual ~psi_initializer();

        // shared methods
        /// @brief allocate memory for psi
        /// @param only_psig if true, only allocate memory for psig, not psi
        /// @return pointer to psi, memory allocated
        /// @note whether use template for psig or not, it is std::complex<double> because psi is.
        psi::Psi<std::complex<double>>* allocate(bool only_psig = false);

        /// @brief get method of initializing psi
        /// @return the method
        std::string get_method() const { return this->method; }
        /// @brief get complement number of bands
        /// @return nbands_complem
        int get_nbands_complem() const { return this->nbands_complem; }

        /// @brief set method manually
        /// @param method_in initialization method
        void set_method(std::string method_in) { this->method = method_in; }
        /// @brief set number of complementary bands
        /// @param nbands_in nbands_complem
        void set_nbands_complem(int nbands_in) { this->nbands_complem = nbands_in; }

        // random to complement bands not initialized by pswfc or nao, therefore it is a basic function, or psi_initializer_random will be inherented by all other methods.
        /// @brief kernel to generate and assign random number for psi
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        void random_t(T* psi, const int iw_start, const int iw_end, const int ik);

        // random
        /// @brief wrapper of random_t
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        virtual void random(T* psi, const int iw_start, const int iw_end,
                            const int ik) { ModuleBase::WARNING_QUIT("psi_initializer::random", "Polymorphism error"); }
        #ifdef __MPI
        /// @brief (about planewaves distribution) from stick mapping to pool
        /// @param stick 
        /// @param ir 
        /// @param out 
        void stick_to_pool(Real* stick, const int& ir, Real* out) const;
        #endif

        // mutual methods, virtual, will be implemented differently in derived classes
        /// @brief create table for storing calculated overlap between pseudowavefunction/numerical orbitals with spherical bessel function
        virtual void create_ovlp_Xjlq() { ModuleBase::WARNING_QUIT("psi_initializer::create_ovlp_table", "Polymorphism error"); }
        /// @brief calculate psi in planewave representation
        /// @param psi psi
        /// @param ik index of kpoint
        virtual psi::Psi<T, Device>* cal_psig(int ik) = 0;

        /// @brief for variables can be only initialized for once.
        /// @param p_pspot_nl_in (for atomic) interfaces to pseudopot_cell_vnl object, in GlobalC now
        /// @attention if one variable is necessary for all methods, initialize it in constructor, not here.
        virtual void initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in = nullptr) = 0;
        // atomic
        /// @brief setter of pseudopotential files, useful when init_wfc = atomic
        virtual void set_pseudopot_files(std::string* pseudopot_files) { ModuleBase::WARNING_QUIT("psi_initializer::set_pseudopot_files", "Polymorphism error"); }
        /// @brief normalize pseudo wavefunction
        /// @param n_rgrid level of realspace grid
        /// @param pswfc pseudowavefunction read from pseudopotential file
        /// @param rgrid realspace grid read from pseudopotential file
        virtual void normalize_pswfc(int n_rgrid, double* pswfc, double* rgrid) { ModuleBase::WARNING_QUIT("psi_initializer::normalize_pswfc", "Polymorphism error"); }
        /// @brief calculate cos(arg)+isin(arg)
        /// @param arg argument
        /// @param mode if 1, return cos(arg), 0, return cos(arg)+isin(arg), -1, return sin(arg)
        /// @return it depends
        virtual std::complex<double> phase_factor(double arg, int mode = 0) { ModuleBase::WARNING_QUIT("psi_initializer::phase_factor", "Polymorphism error"); return std::complex<double>(0.0,0.0);}

        /// @brief calculate overlap table between pseudowavefunction and spherical bessel function
        virtual void cal_ovlp_pswfcjlq() { ModuleBase::WARNING_QUIT("psi_initializer::calc_ovlp_pswfcjlq", "Polymorphism error"); }
        // nao
        /// @brief setter of numerical orbital files, useful when init_wfc = nao
        virtual void set_orbital_files(std::string* orbital_files) { ModuleBase::WARNING_QUIT("psi_initializer::set_orbital_files", "Polymorphism error"); }
        /// @brief calculate overlap between numerical orbital and spherical bessel function
        virtual void cal_ovlp_flzjlq() { ModuleBase::WARNING_QUIT("psi_initializer::cal_ovlp_flzjlq", "Polymorphism error"); }
        // atomic+random
        // nao+random
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

        /// @brief setter of p_ucell
        /// @param p_ucell_in UnitCell pointer
        void set_interface_ucell(UnitCell* p_ucell_in) { this->p_ucell = p_ucell_in; }
        /// @brief getter of p_ucell
        /// @return this->p_ucell
        UnitCell* get_interface_ucell() const { return this->p_ucell; }
        #ifdef __MPI
        /// @brief setter of p_parakpts
        /// @param p_parakpts_in Parallel_Kpoints pointer
        void set_interface_parakpts(Parallel_Kpoints* p_parakpts_in) { this->p_parakpts = p_parakpts_in; }
        /// @brief getter of p_parakpts
        /// @return this->p_parakpts
        Parallel_Kpoints* get_interface_parakpts() const { return this->p_parakpts; }
        #endif
        /// @brief setter of p_pspot_nl
        /// @param p_pspot_nl_in pseudopot_cell_vnl pointer
        void set_interface_pspot_nl(pseudopot_cell_vnl* p_pspot_nl_in) { this->p_pspot_nl = p_pspot_nl_in; }
        /// @brief getter of p_pspot_nl
        /// @return this->p_pspot_nl
        pseudopot_cell_vnl* get_interface_pspot_nl() const { return this->p_pspot_nl; }
        /// @brief setter of sf
        /// @param sf_in Structure_Factor pointer
        void set_interface_sf(Structure_Factor* sf_in) { this->sf = sf_in; }
        /// @brief getter of sf
        /// @return this->sf
        Structure_Factor* get_interface_sf() const { return this->sf; }
        /// @brief setter of pw_wfc
        /// @param pw_wfc_in ModulePW::PW_Basis_K pointer
        void set_interface_pw_wfc(ModulePW::PW_Basis_K* pw_wfc_in) { this->pw_wfc = pw_wfc_in; }
        /// @brief getter of pw_wfc
        /// @return this->pw_wfc
        ModulePW::PW_Basis_K* get_interface_pw_wfc() const { return this->pw_wfc; }
        /// @brief setter of random_seed
        /// @param random_seed_in new value of random_seed
        void set_random_seed(const int random_seed_in) { this->random_seed = random_seed_in; }
        /// @brief getter of random_seed
        /// @return this->random_seed
        int get_random_seed() const { return this->random_seed; }
        /// @brief setter of mem_saver
        /// @param mem_saver_in new value of mem_saver
        void set_mem_saver(const int mem_saver_in) { this->mem_saver = mem_saver_in; }
        /// @brief getter of mem_saver
        /// @return this->mem_saver
        int get_mem_saver() const { return this->mem_saver; }
        /// @brief setter of initialized
        /// @param initialized_in new value of initialized
        void set_initialized(bool initialized_in) { this->initialized = initialized_in; }
        /// @brief getter of initialized
        /// @return this->initialized
        bool get_initialized() const { return this->initialized; }
        // member variables
        /// @brief interface to the psi::Psi data structure class
        psi::Psi<T, Device>* psig = nullptr;

        /// @brief cast from std::complex<double> to float
        /// @tparam U float placeholder
        /// @param in psi value to cast
        /// @return float psi value
        template <typename U>
        typename std::enable_if<std::is_same<U, float>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return static_cast<float>(in.real());
        }
        /// @brief cast from std::complex<double> to double
        /// @tparam U double placeholder
        /// @param in psi value to cast
        /// @return double psi value
        template <typename U>
        typename std::enable_if<std::is_same<U, double>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return static_cast<double>(in.real());
        }
        /// @brief cast from std::complex<double> to std::complex<float>
        /// @tparam U std::complex<float> placeholder
        /// @param in psi value to cast
        /// @return std::complex<float> psi value
        template <typename U>
        typename std::enable_if<std::is_same<U, std::complex<float>>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return std::complex<float>(static_cast<float>(in.real()), static_cast<float>(in.imag()));
        }
        /// @brief cast from std::complex<double> to std::complex<double>
        /// @tparam U std::complex<double> placeholder
        /// @param in psi value to cast
        /// @return std::complex<double> psi value
        template <typename U>
        typename std::enable_if<std::is_same<U, std::complex<double>>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return std::complex<double>(in.real(), in.imag());
        }
        Real norm2(const std::complex<Real> in)
        {
            return in.real()*in.real() + in.imag()*in.imag();
        }
        Real norm2(const Real in)
        {
            return in*in;
        }
        /*
        template <typename U>
        void cast_right_to_left(U& left, const std::complex<double> right)
        {
            if(std::is_same<U, float>::value) *(*float)left = static_cast<float>(right.real());
            else if(std::is_same<U, double>::value) *(*double)left = static_cast<double>(right.real());
            else if(std::is_same<U, std::complex<float>>::value) *(*std::complex<float>)left = std::complex<float>(static_cast<float>(right.real()), static_cast<float>(right.imag()));
            else if(std::is_same<U, std::complex<double>>::value) *(*std::complex<double>)left = std::complex<double>(right.real(), right.imag());
            else ModuleBase::WARNING_QUIT("psi_initializer::cast_right_to_left", "type error");
        }
        */
    protected:
        // interfaces
        // ATTENTION: DO NOT USE DELETE ON THESE POINTERS
        // normal interfaces
        /// @brief interface to the Structure_Factor method class
        Structure_Factor* sf = nullptr;
        /// @brief interface to the PW_Basis_K data structure class
        ModulePW::PW_Basis_K* pw_wfc = nullptr;
        // interfaces designed to get rid of dependence troubles of GlobalC in unittest
        /// @brief interface to the UnitCell data carrier class, used in all methods
        UnitCell* p_ucell = nullptr;
        #ifdef __MPI
        /// @brief interface to the Parallel_Kpoints method class, used in all methods
        Parallel_Kpoints* p_parakpts = nullptr;
        #endif
        /// @brief interface to the pseudopot_cell_vnl data carrier class, used in atomic
        pseudopot_cell_vnl* p_pspot_nl = nullptr;

        int random_seed = 1; // random seed
        // tool interfaces
        /// @brief method of Spherical Bessel Transformation
        ModuleBase::SphericalBesselTransformer sbt; // useful for atomic-like methods
    private:
        // basic properties
        int mem_saver = 0; // will deprecated this variable soon
        std::string method = "none";
        // non-random case
        int nbands_complem = 0;
        // random
        int* ixy2is;

        bool initialized = false; // whether initialized or not
        // atomic+random or nao+random
        double random_mix = 0;
};
#endif