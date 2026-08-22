#ifndef PSI_INIT_ATOM_RAND_H
#define PSI_INIT_ATOM_RAND_H
#include "psi_init_atomic.h"

/*
Psi (planewave based wavefunction) initializer: atomic+random
*/
template <typename T>
class psi_init_atomic_random : public psi_init_atomic<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_init_atomic_random()
    {
        this->method_ = "atomic+random";
        this->mixing_coef_ = 0.05;
    }
    ~psi_init_atomic_random(){};

    virtual void initialize(const Structure_Factor* sf,
                            const ModulePW::PW_Basis_K* pw_wfc,
                            const UnitCell* p_ucell,
                            const std::vector<int>& ik2iktot,
                            const int& random_seed,
                            const int& rank,
                            const int& npol,
                            const int& nbands) override;

    virtual void init_psig(T* psig, const int& ik) override;

  private:
};
#endif