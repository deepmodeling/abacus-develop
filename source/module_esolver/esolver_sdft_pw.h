#ifndef ESOLVER_SDFT_PW_H
#define ESOLVER_SDFT_PW_H

#include "esolver_ks_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
#include "module_hamilt_pw/hamilt_stodft/sto_hchi.h"

namespace ModuleESolver
{

class ESolver_SDFT_PW: public ESolver_KS_PW<std::complex<double>>
{
public:
    ESolver_SDFT_PW();
    ~ESolver_SDFT_PW();
    void Init(Input &inp, UnitCell &cell) override;
    double cal_Energy() override;
    void cal_Force(ModuleBase::matrix& force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;
public:
    Stochastic_WF stowf;
    Stochastic_Iter stoiter;

protected:
    virtual void beforescf(const int istep) override; 
    // virtual void eachiterinit(int iter) override; 
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void nscf() override;
    virtual void othercalculation(const int istep) override;
    virtual void eachiterfinish(const int iter) override; 
    virtual void afterscf(const int istep) override;
    virtual void postprocess() override;

public:
    //Stochastic Kubo-Greenwood
    void sKG(const int nche_KG, const double fwhmin, const double wcut, 
             const double dw_in, const double dt_in, const int nbatch);
    //calculate DOS
    void caldos(const int nche_dos, const double sigmain, 
            const double emin, const double emax, const double de, const int npart);

private:
    int nche_sto; //norder of Chebyshev
    void check_che(const int nche_in);


};

}//namespace ModuleESolver

//temporary setting: removed GlobalC but not breaking design philosophy
namespace GlobalTemp
{

extern const ModuleBase::matrix* veff;

}

#endif