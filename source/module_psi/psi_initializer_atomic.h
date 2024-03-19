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
        psi_initializer_atomic() {this->set_method("atomic");}
        ~psi_initializer_atomic() {};

        psi::Psi<T, Device>* cal_psig(int ik) override;

        #ifdef __MPI // MPI additional implementation
        /// @brief initialize the psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,                      //< structure factor
                                ModulePW::PW_Basis_K*,                  //< planewave basis
                                UnitCell*,                              //< unit cell
                                Parallel_Kpoints*,                      //< parallel kpoints
                                const int& = 1,                         //< random seed
                                pseudopot_cell_vnl* = nullptr,          //< nonlocal pseudopotential
                                const int& = 0) override;               //< MPI rank
        #else
        /// @brief serial version of initialize function, link psi_initializer with external data and methods
        virtual void initialize(Structure_Factor*,                      //< structure factor
                                ModulePW::PW_Basis_K*,                  //< planewave basis
                                UnitCell*,                              //< unit cell
                                const int& = 1,                         //< random seed
                                pseudopot_cell_vnl* = nullptr) override;//< nonlocal pseudopotential
        #endif

        void allocate_table() override;
        void tabulate() override;

        std::vector<std::string> pseudopot_files() const { return pseudopot_files_; }

    private:
        std::vector<std::string> pseudopot_files_;
        ModuleBase::realArray ovlp_pswfcjlq_;
};
#endif