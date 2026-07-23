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

    /// @brief initialize the psi_init with external data and methods
    virtual void initialize(const Structure_Factor*,             //< structure factor
                            const ModulePW::PW_Basis_K*,         //< planewave basis
                            const UnitCell*,                     //< unit cell
                            const std::vector<int>& = {},        //< ik2iktot: local->global k-point mapping
                            const int& = 0,                      //< nkstot: total number of k-points
                            const int& = 1,                      //< random seed
                            const int& = 0,                      //< lmaxkb: max angular momentum for non-local projectors
                            const int& = 0,                      //< MPI rank
                            const int& = 1,                      //< npol
                            const int& = 1) override;            //< nbands

    virtual void init_psig(T* psig, const int& ik) override;
};
#endif