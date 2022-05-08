#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"


namespace ModuleESolver
{

class ESolver_KS_PW: public ESolver_KS
{
public:
    ESolver_KS_PW();
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;
    void postprocess() override;

protected:
    virtual void beforescf(const int istep) override; 
    virtual void eachiterinit(const int istep, const int iter) override; 
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void updatepot(const int istep, const int iter, const bool conv) override;
    virtual void eachiterfinish(const int iter, const bool conv) override; 
    virtual void afterscf(const int iter, const bool conv) override;

    // <Temporary> Get wavefunctions and eigen energies. 
    // It should be replaced by diag class in HSolver module in the future
    void c_bands(const int istep, const int iter);

    // It copies the function in Threshold_Elec class.
    // After all ESolver, HSolver are constructed, Class Electrons and Threshold_Elec should be deleted.
    void print_eigenvalue(std::ofstream &ofs);


    // ESolver_KS_PW(bool use_sdft)
    // {
    //     if(use_sdft)    
    //     {
    //         this->ph2e=new H2E_SDFT();
    //     }
    //     else           
    //     {
    //         this->ph2e=new H2E_PW();
    //     }
    //     this->pes= new Estate_PW();
    //     this->phamilt=new Hamilt_PW();
    // }
    
    // Basis_PW basis_pw;
    // Init(Inputs &inp, Cell &cel)
    // {
        
    //     basis_pw.init(inp, cel);

    //     pes->init(inp, cel, basis_pw); 
        
    //     phamilt->init(bas); 
    //     phamilt->initpot(cel, pes); 
        
    //     ph2e->init(h, pes); 
    // }
};
}
#endif