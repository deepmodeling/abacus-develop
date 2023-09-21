#ifndef PSI_INITIALIZER_ATOMIC_H
#define PSI_INITIALIZER_ATOMIC_H
#include "psi_initializer.h"
#include "module_base/realarray.h"

/*
Psi (planewave based wavefunction) initializer: atomic
*/
class psi_initializer_atomic : public psi_initializer
{
    public:

        /// @brief constructor of psi initializer (atomic)
        /// @param sf_in Structure factor interface, link ESolver
        /// @param pw_wfc_in ModulePW::PW_Basis_K interface, link ESolver
        psi_initializer_atomic(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        /// @brief default destructor
        ~psi_initializer_atomic();

        // methods
        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<std::complex<double>>* cal_psig(int ik) override;

        // setters

        /// @brief setter of pseudpotential filenames
        /// @param pseudopot_files pseudpotential filenames organized in an array
        void set_pseudopot_files(std::string* pseudopot_files);
        // I wont write a function to set ovlp_pswfcjlq, it is totally useless

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
        /// @brief calculate the overlap between pseudo atomic wavefunctions and planewave basis
        void cal_ovlp_pswfcjlq();

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
#endif