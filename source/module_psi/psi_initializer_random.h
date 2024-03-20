#ifndef PSI_INITIALIZER_RANDOM_H
#define PSI_INITIALIZER_RANDOM_H

#include "psi_initializer.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"

/*
Psi (planewave based wavefunction) initializer: random method
*/
template <typename T, typename Device>
class psi_initializer_random : public psi_initializer<T, Device>
{
    private:
        using Real = typename GetTypeReal<T>::type;
    public:
        psi_initializer_random() {this->set_method("random");};
        ~psi_initializer_random() {};
        /// @brief write random number to psi in certain range specified by ik, iw_start, iw_end
        void random(T* psi,                 //< psi
                    const int iw_start,     //< iw_start, starting band index of present kpoint
                    const int iw_end,       //< iw_end, ending band index of present kpoint
                    const int ik) override; //< ik, kpoint index
        /// @brief calculate and output planewave wavefunction
        /// @param ik kpoint index
        /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
        psi::Psi<T, Device>* cal_psig(int ik) override;
        #ifdef __MPI // MPI additional implementation
        /// @brief initialize the psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,              //< structure factor
                                ModulePW::PW_Basis_K*,          //< planewave basis
                                UnitCell*,                      //< unit cell
                                Parallel_Kpoints*,              //< parallel kpoints
                                const int& = 1,                 //< random seed
                                pseudopot_cell_vnl* = nullptr,  //< nonlocal pseudopotential
                                const int& = 0) override;       //< MPI rank
        #else
        /// @brief serial version of initialize function, link psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,                          //< structure factor
                                ModulePW::PW_Basis_K*,                      //< planewave basis
                                UnitCell*,                                  //< unit cell
                                const int& = 1,                             //< random seed
                                pseudopot_cell_vnl* = nullptr) override   ; //< nonlocal pseudopotential
        #endif
        virtual void tabulate() override {};
};
#endif