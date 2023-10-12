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
//#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
//#endif
template<typename T, typename Device>
class RepIn
{
    public:
        RepIn() = delete;
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
        virtual ~RepIn() { };

        virtual void cal_psig(const psi::Psi<T, Device>* psig) = 0;
        std::string get_representation() const { return this->representation; };
        void set_kpoint(int ik_in) { this->ik = ik_in; }
        // methods
        // mutual methods, virtual, will be implemented differently in derived classes
        /// @brief create table for storing calculated overlap between pseudowavefunction/numerical orbitals with spherical bessel function
        virtual void create_ovlp_Xjlq() { ModuleBase::WARNING_QUIT("psi_initializer::create_ovlp_table", "Polymorphism error"); }

        // pao
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

        /// @brief setter of p_ucell
        /// @param p_ucell_in UnitCell pointer
        void set_interface_ucell(UnitCell* p_ucell_in) { this->p_ucell = p_ucell_in; }
        /// @brief getter of p_ucell
        /// @return this->p_ucell
        UnitCell* get_interface_ucell() const { return this->p_ucell; }
        //#ifdef __MPI
        /// @brief setter of p_parakpts
        /// @param p_parakpts_in Parallel_Kpoints pointer
        void set_interface_parakpts(Parallel_Kpoints* p_parakpts_in) { this->p_parakpts = p_parakpts_in; }
        /// @brief getter of p_parakpts
        /// @return this->p_parakpts
        Parallel_Kpoints* get_interface_parakpts() const { return this->p_parakpts; }
        //#endif
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
    protected:
        int ik = 0; // kpoint index
        // interfaces
        Structure_Factor* p_sf;
        ModulePW::PW_Basis_K* pw_wfc;
        UnitCell* p_ucell;
        Parallel_Kpoints* p_parakpts;
        pseudopot_cell_vnl* p_pspot_nl;
        // numerical algorithm support
        ModuleBase::SphericalBesselTransformer sbt;
};

#endif // REPIN_H