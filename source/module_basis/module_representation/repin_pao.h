#ifndef REPIN_PAO_H
#define REPIN_PAO_H

#include "module_basis/module_representation/repin.h"

template<typename T, typename Device>
class RepIn_PAO : public RepIn<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        RepIn_PAO() = delete;
        #ifdef __MPI
        /// @brief constructor of RepIn
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_parakpts_in link to Parallel_Kpoints
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn_PAO(Structure_Factor* sf_in, 
                  ModulePW::PW_Basis_K* pw_wfc_in, 
                  UnitCell* p_ucell_in, 
                  Parallel_Kpoints* p_parakpts_in, 
                  pseudopot_cell_vnl* p_pspot_nl_in);
        #else
        /// @brief constructor of RepIn
        /// @param sf_in link to Structure_Factor
        /// @param pw_wfc_in link to ModulePW::PW_Basis_K
        /// @param p_ucell_in link to UnitCell
        /// @param p_pspot_nl_in link to pseudopot_cell_vnl
        RepIn_PAO(Structure_Factor* sf_in, 
                  ModulePW::PW_Basis_K* pw_wfc_in, 
                  UnitCell* p_ucell_in, 
                  pseudopot_cell_vnl* p_pspot_nl_in);
        #endif
        ~RepIn_PAO();
        /*
            central function
        */
        /// @brief calculate pw representation of basis function of input wavefunction
        /// @param psig basis function in pw representation
        void cal_psig(psi::Psi<T, Device>* psig) override;
        /*
            initialization and subroutines
        */
        /// @brief call create_ovlp_pswfcjlq, cal_ovlp_pswfcjlq
        /// @param pseudopot_files array storing pseudopotential files, not used in PAO
        void initialize(std::string* pseudopot_files) override;
        /// @brief allocate memory for ovlp_pswfcjlq and initialize all elements to 0
        void create_ovlp_pswfcjlq();
        /// @brief calculate the overlap between pseudo atomic wavefunctions and planewave basis
        void cal_ovlp_pswfcjlq();
        /*
            elemental operations
        */
        /// @brief specialized normalization of wfc function
        /// @param n_rgrid number of grid points in realspace
        /// @param pswfc pseudo wavefunction in pseudopotential files
        /// @param rgrid realspace grid points, r1, r2, ...
        void normalize_pswfc(int n_rgrid, double* pswfc, double* rgrid);
        /// @brief simple unitary phase factor
        /// @param arg the argument of the phase factor
        /// @param mode +1 for real part, -1 for imaginary part, 0 for the whole
        /// @return the phase factor
        std::complex<double> phase_factor(double arg, int mode = 0);
        /*
            getters
        */
        /// @brief getter of pseudpotential files list
        /// @return pseudopotential files list
        std::vector<std::string> get_pseudopot_files() const { return pseudopot_files; }
        /// @brief getter of matrix of overlap between pseudo wavefunction and spherical bessel function
        /// @return ovlp_pswfcjlq
        ModuleBase::realArray get_ovlp_pswfcjlq() const { return ovlp_pswfcjlq; }
    private:
        std::vector<std::string> pseudopot_files;
        ModuleBase::realArray ovlp_pswfcjlq;
};

#endif // REPIN_PAO_H