#ifndef PSI_INITIALIZER_H
#define PSI_INITIALIZER_H
// data structure support
#include "module_psi/psi.h" // for psi data structure
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "module_basis/module_pw/pw_basis_k.h" // for kpoint related data structure
// numerical algorithm support
#include "module_base/spherical_bessel_transformer.h" // for spherical bessel transform
#ifdef __MPI
// MPI support
#include <mpi.h>
#include "module_base/parallel_global.h"
#include "module_cell/parallel_kpoints.h"
#endif
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
class psi_initializer
{
    public:
        /// @brief default constructor of psi initializer
        psi_initializer() : sf(nullptr), pw_wfc(nullptr) { };
        /// @brief parameterized constructor of psi initializer
        /// @param sf_in structure factor pointer
        /// @param pw_wfc_in pw_basis_k pointer
        psi_initializer(Structure_Factor* sf_in, ModulePW::PW_Basis_K* pw_wfc_in);
        /// @brief destructor
        virtual ~psi_initializer();

        /// @brief allocate memory for psi
        /// @return pointer to psi, memory allocated
        psi::Psi<std::complex<double>>* allocate();

        /// @brief calculate psi in planewave representation
        /// @param psi psi
        /// @param ik index of kpoint
        virtual psi::Psi<std::complex<double>>* cal_psig(int ik) = 0;

        /// @brief get method of initializing psi
        /// @return the method
        std::string get_method() const { return this->method; }
        /// @brief get complement number of bands
        /// @return nbands_complem
        int get_nbands_complem() const { return this->nbands_complem; }

        void print_status(psi::Psi<std::complex<double>>& psi) const;

        /// @brief set method manually
        /// @param method_in initialization method
        void set_method(std::string method_in) { this->method = method_in; }
        /// @brief set number of complementary bands
        /// @param nbands_in nbands_complem
        void set_nbands_complem(int nbands_in) { this->nbands_complem = nbands_in; }
        /// @brief output psig to file, given number of kpoints, bands and basis, diagonalization method. Note: because it is complex number, therefore every number has format (real, imag)
        void write_psig() const;
        /// @brief output psig at given kpoint to file
        void write_psig(int ik) const;
        // virtual functions, will be implemented in derived classes

        // random to complement bands not initialized by pswfc or nao, therefore it is a basic function, or psi_initializer_random will be inherented by all other methods.
        /// @brief kernel to generate and assign random number for psi
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        /// @param wfc_basis because do not want to change contents in it, this parameter is reserved. should give this->pw_wfc
        void random_t(std::complex<double>* psi, const int iw_start, const int iw_end, const int ik, const ModulePW::PW_Basis_K* wfc_basis);
        
        // random
        /// @brief wrapper of random_t
        /// @param psi for psi::Psi<FPTYPE, Device>, first use psi.fix(ik), then psi.get_pointer() to pass to parameter list
        /// @param iw_start index of wavefunction (start), the index of first band to initialize
        /// @param iw_end index of wavefunction (end)
        /// @param ik index of kpoint
        /// @param wfc_basis because do not want to change contents in it, this parameter is reserved. should give this->pw_wfc
        virtual void random(std::complex<double>* psi, const int iw_start, const int iw_end,
                            const int ik, const ModulePW::PW_Basis_K* wfc_basis) { ModuleBase::WARNING_QUIT("psi_initializer::random", "Polymorphism error"); }
        #ifdef __MPI
        /// @brief (about planewaves distribution) from stick mapping to pool
        /// @param stick 
        /// @param ir 
        /// @param out 
        /// @param wfc_basis because do not want to change contents in it, this parameter is reserved. should give this->pw_wfc
        void stick_to_pool(double* stick, const int& ir, double* out, const ModulePW::PW_Basis_K* wfc_basis) const;
        #endif
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
        // member variables
        /// @brief interface to the psi::Psi data structure class
        psi::Psi<std::complex<double>>* psig;
        /// @brief interface to the Structure_Factor method class
        Structure_Factor* sf;
        /// @brief interface to the PW_Basis_K data structure class
        ModulePW::PW_Basis_K* pw_wfc;
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

        // atomic+random or nao+random
        double random_mix = 0;
};
#endif