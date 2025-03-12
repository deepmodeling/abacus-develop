#ifndef PSI_INITIALIZER_ATOMIC_RANDOM_H
#define PSI_INITIALIZER_ATOMIC_RANDOM_H
#include "module_cell/parallel_kpoints.h"
#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "psi_initializer_atomic.h"

/*
Psi (planewave based wavefunction) initializer: atomic+random
*/
template <typename T>
class psi_initializer_atomic_random : public psi_initializer_atomic<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_initializer_atomic_random()
    {
        this->method_ = "atomic+random";
        this->mixing_coef_ = 0.05;
    }
    ~psi_initializer_atomic_random(){};

    /// @brief initialize the psi_initializer with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const K_Vectors*,                    //< kpoints
                            const int& = 1,                      //< random seed
                            const pseudopot_cell_vnl* = nullptr, //< nonlocal pseudopotential
                            const int& = 0) override;            //< MPI rank

    virtual void init_psig(T* psig, const int& ik) override;

  private:
};
#endif