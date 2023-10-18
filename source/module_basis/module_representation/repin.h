#ifndef REPIN_H
#define REPIN_H
// data structure support
#include <string>
#include "module_psi/psi.h"
// numerical algorithm support
#include "module_base/spherical_bessel_transformer.h" // for spherical bessel transform
// for RepIn and RepOut initialization
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#endif
// template programming
#include "module_base/macros.h"
/*
    RepIn: A class transforming basis function to pw representation

    Usage: cal_psig(psig);
*/

template<typename T, typename Device>
class RepIn
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
#ifdef __MPI
        /// @brief constructor of RepIn
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_parakpts_in link to Parallel_Kpoints
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn(Structure_Factor* sf_in, 
              ModulePW::PW_Basis_K* pw_wfc_in, 
              UnitCell* p_ucell_in, 
              Parallel_Kpoints* p_parakpts_in, 
              pseudopot_cell_vnl* p_pspot_nl_in): 
              p_sf(sf_in), 
              pw_wfc(pw_wfc_in), 
              p_ucell(p_ucell_in), 
              p_parakpts(p_parakpts_in),
              p_pspot_nl(p_pspot_nl_in) { };
#else
        /// @brief constructor of RepIn
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn(Structure_Factor* sf_in, 
              ModulePW::PW_Basis_K* pw_wfc_in, 
              UnitCell* p_ucell_in, 
              pseudopot_cell_vnl* p_pspot_nl_in): 
              p_sf(sf_in), 
              pw_wfc(pw_wfc_in), 
              p_ucell(p_ucell_in), 
              p_pspot_nl(p_pspot_nl_in) { };
#endif
        virtual ~RepIn() { };

        /* ------------------
            central function
           ------------------ */
        /// @brief do transform on basis function from one representation to pw
        /// @param psig pw representation of basis function of representation
        virtual void cal_psig(psi::Psi<T, Device>* psig) = 0;
        /* ------------------
            initialization (and subroutines)
           ------------------ */
        /// @brief initialization of pw representation, will allocate overlap table, read files, etc.
        virtual void initialize(std::string* filenames = nullptr) { ModuleBase::WARNING_QUIT("RepIn::initialization", "Polymorphism error"); }
        /// @brief getter of representation input
        /// @return representation_in
        std::string get_representation() const { return this->representation; };
        /// @brief setter of kpoint index
        /// @param ik_in kpoint index
        void set_kpoint(int ik_in) { this->ik = ik_in; }
        /* ------------------
            interfaces setting
           ------------------ */
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
        void set_interface_sf(Structure_Factor* sf_in) { this->p_sf = sf_in; }
        /// @brief getter of sf
        /// @return this->sf
        Structure_Factor* get_interface_sf() const { return this->p_sf; }
        /// @brief setter of pw_wfc
        /// @param pw_wfc_in ModulePW::PW_Basis_K pointer
        void set_interface_pw_wfc(ModulePW::PW_Basis_K* pw_wfc_in) { this->pw_wfc = pw_wfc_in; }
        /// @brief getter of pw_wfc
        /// @return this->pw_wfc
        ModulePW::PW_Basis_K* get_interface_pw_wfc() const { return this->pw_wfc; }
        /* ------------------
            multiple float type and device support
           ------------------ */
        /// @brief cast from std::complex<double> to float
        /// @tparam U float placeholder
        /// @param in psi value to cast
        /// @return float psi value
        template <typename U>
        static typename std::enable_if<std::is_same<U, float>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return static_cast<float>(in.real());
        }
        /// @brief cast from std::complex<double> to double
        /// @tparam U double placeholder
        /// @param in psi value to cast
        /// @return double psi value
        template <typename U>
        static typename std::enable_if<std::is_same<U, double>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return static_cast<double>(in.real());
        }
        /// @brief cast from std::complex<double> to std::complex<float>
        /// @tparam U std::complex<float> placeholder
        /// @param in psi value to cast
        /// @return std::complex<float> psi value
        template <typename U>
        static typename std::enable_if<std::is_same<U, std::complex<float>>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return std::complex<float>(static_cast<float>(in.real()), static_cast<float>(in.imag()));
        }
        /// @brief cast from std::complex<double> to std::complex<double>
        /// @tparam U std::complex<double> placeholder
        /// @param in psi value to cast
        /// @return std::complex<double> psi value
        template <typename U>
        static typename std::enable_if<std::is_same<U, std::complex<double>>::value, U>::type cast_to_T(const std::complex<double> in)
        {
            return std::complex<double>(in.real(), in.imag());
        }
    protected:
        int ik = 0; // kpoint index
        // interfaces
        Structure_Factor* p_sf = nullptr;
        ModulePW::PW_Basis_K* pw_wfc = nullptr;
        UnitCell* p_ucell = nullptr;
#ifdef __MPI
        Parallel_Kpoints* p_parakpts = nullptr;
#endif
        pseudopot_cell_vnl* p_pspot_nl = nullptr;
        // numerical algorithm support
        ModuleBase::SphericalBesselTransformer sbt;

        std::string representation = "unknown"; // representation of input psi
};

#endif // REPIN_H