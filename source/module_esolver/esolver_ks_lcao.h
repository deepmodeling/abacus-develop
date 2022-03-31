#ifndef ESOLVER_KS_LCAO_H
#define ESOLVER_KS_LCAO_H
#include "./esolver_ks.h"

#include "../src_lcao/LOOP_elec.h"

namespace ModuleESolver
{

class ESolver_KS_LCAO: public ESolver_KS
{
public:
    ESolver_KS_LCAO()
    {
        tag = "ESolver_KS_LCAO";
    }
    void Init(Input& inp, UnitCell_pseudo& cell) override;
    
    void Run(int istep,
        Record_adj& ra,
        Local_Orbital_Charge& loc,
        Local_Orbital_wfc& lowf,
        LCAO_Hamilt& uhm) override;
    void Run(int istep, UnitCell_pseudo& cell) override {};
    
    void cal_Energy(energy& en) override;
    void cal_Force(ModuleBase::matrix &force) override;
    void cal_Stress(ModuleBase::matrix& stress) override;

private:
    LOOP_elec LOE;
};

///Basis_lcao
///ORB_control orb_con;


}
#endif