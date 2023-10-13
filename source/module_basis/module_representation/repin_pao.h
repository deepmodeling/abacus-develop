#ifndef REPIN_PAO_H
#define REPIN_PAO_H

#include "module_basis/module_representation/repin.h"

template<typename T, typename Device>
class RepIn_PAO : public RepIn<T, Device>
{
    public:
        RepIn_PAO() = delete;
        RepIn_PAO(Structure_Factor* sf_in, 
                  ModulePW::PW_Basis_K* pw_wfc_in, 
                  UnitCell* p_ucell_in, 
                  Parallel_Kpoints* p_parakpts_in, 
                  pseudopot_cell_vnl* p_pspot_nl_in) = default;
        ~RepIn_PAO();
        void cal_psig(const psi::Psi<T, Device>* psig) override;

        /// @brief allocate memory for ovlp_pswfcjlq and initialize all elements to 0
        void create_ovlp_Xjlq() override;
        /// @brief specialized normalization of wfc function
        /// @param n_rgrid number of grid points in realspace
        /// @param pswfc pseudo wavefunction in pseudopotential files
        /// @param rgrid realspace grid points, r1, r2, ...
        void normalize_pswfc(int n_rgrid, double* pswfc, double* rgrid) override;
        /// @brief simple unitary phase factor
        /// @param arg the argument of the phase factor
        /// @param mode +1 for real part, -1 for imaginary part, 0 for the whole
        /// @return the phase factor
        std::complex<double> phase_factor(double arg, int mode = 0) override;
        /// @brief calculate the overlap between pseudo atomic wavefunctions and planewave basis
        void cal_ovlp_pswfcjlq() override;

        void representation_init(const std::string* pseudopot_files) override;
        // historically left functions
        // getters

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