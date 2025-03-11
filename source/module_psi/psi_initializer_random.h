#ifndef PSI_INITIALIZER_RANDOM_H
#define PSI_INITIALIZER_RANDOM_H

#include "module_hamilt_pw/hamilt_pwdft/VNL_in_pw.h"
#include "psi_initializer.h"

/*
Psi (planewave based wavefunction) initializer: random method
*/
template <typename T>
class psi_initializer_random : public psi_initializer<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_initializer_random()
    {
        this->method_ = "random";
    };
    ~psi_initializer_random(){};
    /// @brief calculate and output planewave wavefunction
    /// @param ik kpoint index
    /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
    virtual void init_psig(T* psig, const int& ik) override;
    /// @brief initialize the psi_initializer with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const K_Vectors*,                    //< kpoints
                            const int& = 1,                      //< random seed
                            const pseudopot_cell_vnl* = nullptr, //< nonlocal pseudopotential
                            const int& = 0) override;            //< MPI rank
};
#endif