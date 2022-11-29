#ifndef H_TDDFT_PW_H
#define H_TDDFT_PW_H

#include "pot_base.h"

namespace elecstate
{

class H_TDDFT_pw : public PotBase
{
    public:
    H_TDDFT_pw(
        const ModulePW::PW_Basis* rho_basis_in,
        const UnitCell* ucell_in):ucell_(ucell_in)
    {
        this->dynamic_mode = false;
        this->fixed_mode = true;

        this->rho_basis_ = rho_basis_in;
    }
    ~H_TDDFT_pw(){};

    void cal_fixed_v(double* vl_pseudo) override;

    private:
    //internal time-step, 
    //-------hypothesis-------
    //Vext will evolve by time, every time cal_fixed_v() is called, istep++
    //------------------------ 
    static int istep;

    const UnitCell* ucell_ = nullptr;
};

} //namespace elecstate

#endif