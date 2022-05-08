#include "./esolver_ks_pw.h"
namespace ModuleESolver
{

class ESolver_SDFT_PW: public ESolver_KS_PW
{
public:
    ESolver_SDFT_PW();
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix &stress) override;

protected:
    virtual void beforescf(const int istep) override;
    // virtual void eachiterinit(int iter) override; 
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void eachiterfinish(const int iter, const bool conv) override; 
    virtual void afterscf(const int iter, const bool conv) override;
};

}

