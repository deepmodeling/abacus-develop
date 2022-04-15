#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
// #include "Basis_PW.h"
// #include "Estate_PW.h"
// #include "Hamilton_PW.h"
// #include "H2E_pw.h"

#include "../src_pw/electrons.h"

namespace ModuleESolver
{

class ESolver_KS_PW: public ESolver_KS
{
public:
    ESolver_KS_PW()
    {
        tag = "ESolver_KS_PW";
    }
    void Init(Input &inp, UnitCell_pseudo &cell) override;
    void Run(int istep, UnitCell_pseudo& cell) override;
    void Run(int istep,
        Record_adj& ra,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& lowf,
        LCAO_Hamilt& uhm) override {};
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix &stress) override;

//--------------temporary----------------------------
    Electrons elec;
//---------------------------------------------------
    int getiter();

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