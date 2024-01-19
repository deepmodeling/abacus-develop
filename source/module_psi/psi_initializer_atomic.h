#ifndef PSI_INITIALIZER_ATOMIC_H
#define PSI_INITIALIZER_ATOMIC_H
#include "psi_initializer.h"
#include "module_base/realarray.h"

/*
Psi (planewave based wavefunction) initializer: atomic
*/
template <typename T, typename Device>
class psi_initializer_atomic : public psi_initializer<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        #ifdef __MPI
        /// @brief parameterized constructor of psi initializer (with MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param p_parakpts_in interface, link with Parallel_Kpoints GlobalC::Pkpoints
        /// @param random_seed_in random seed
        psi_initializer_atomic(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, Parallel_Kpoints* p_parakpts_in, int random_seed_in = 1);
        #else
        /// @brief parameterized constructor of psi initializer (without MPI support)
        /// @param sf_in interface, link with Structure_Factor ESolver_FP::sf
        /// @param pw_wfc_in interface, link with ModulePW::PW_Basis_K* ESolver_FP::pw_wfc
        /// @param p_ucell_in interface, link with UnitCell GlobalC::ucell
        /// @param random_seed_in random seed
        psi_initializer_atomic(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in, UnitCell* p_ucell_in, int random_seed_in = 1);
        #endif
        /// @brief default destructor
        ~psi_initializer_atomic();

        // methods
        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<T, Device>* cal_psig(int ik) override;

        /// @brief initialize only once, for atomic, it should be, create ovlp_pswfcjlq, calculate ovlp_pswfcjlq
        /// @param p_pspot_nl_in (for atomic) interfaces to pseudopot_cell_vnl object, in GlobalC now
        /// @attention if one variable is necessary for all methods, initialize it in constructor, not here.
        void initialize_only_once(pseudopot_cell_vnl* p_pspot_nl_in) override;
        // setters

        /* I leave this function here for deprecation of UnitCell in the future */
        /// @brief setter of pseudpotential filenames
        /// @param pseudopot_files pseudpotential filenames organized in an array
        //void set_pseudopot_files(std::string* pseudopot_files);
        // I wont write a function to set ovlp_pswfcjlq, it is totally useless

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