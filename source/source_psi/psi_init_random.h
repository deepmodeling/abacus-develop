#ifndef PSI_INIT_RANDOM_H
#define PSI_INIT_RANDOM_H

#include <vector>
#include "psi_base.h"

/*
Psi (planewave based wavefunction) initializer: random method
*/
template <typename T>
class psi_init_random : public psi_base<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_init_random()
    {
        this->method_ = "random";
    };
    ~psi_init_random(){};
    /// @brief calculate and output planewave wavefunction
    /// @param ik kpoint index
    /// @return initialized planewave wavefunction (psi::Psi<std::complex<double>>*)
    virtual void init_psig(T* psig, const int& ik) override;
    virtual void initialize(const Structure_Factor* sf,             //< structure factor
                            const ModulePW::PW_Basis_K* pw_wfc,         //< planewave basis
                            const UnitCell* p_ucell,                     //< unit cell
                            const std::vector<int>& ik2iktot,             //< ik2iktot: local->global k-point mapping
                            const int& random_seed,                      //< random seed
                            const int& rank,                            //< MPI rank
                            const int& npol,                            //< npol
                            const int& nbands) override;                //< nbands
};
#endif