#ifndef PSI_INITIALIZER_NAO_RANDOM_H
#define PSI_INITIALIZER_NAO_RANDOM_H
#include "module_cell/parallel_kpoints.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "psi_initializer_nao.h"

/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital + random method
*/
template <typename T>
class psi_initializer_nao_random : public psi_initializer_nao<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_initializer_nao_random()
    {
        this->method_ = "nao+random";
        this->mixing_coef_ = 0.05;
    };
    ~psi_initializer_nao_random(){};

    /// @brief initialize the psi_initializer with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const Parallel_Kpoints*,             //< parallel kpoints
                            const int& = 1,                //< random seed
                            const pseudopot_cell_vnl* = nullptr, //< nonlocal pseudopotential
                            const int& = 0) override;      //< MPI rank

    virtual void init_psig(T* psig,  const int& ik) override;
};
#endif