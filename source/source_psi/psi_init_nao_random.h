#ifndef PSI_INIT_NAO_RANDOM_H
#define PSI_INIT_NAO_RANDOM_H
#include "psi_init_nao.h"

/*
Psi (planewave based wavefunction) initializer: numerical atomic orbital + random method
*/
template <typename T>
class psi_init_nao_random : public psi_init_nao<T>
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    psi_init_nao_random()
    {
        this->method_ = "nao+random";
        this->mixing_coef_ = 0.05;
    };
    ~psi_init_nao_random(){};

    virtual void init_psig(T* psig, const int& ik) override;
};
#endif