#ifndef CHARGE_EXTRA_H
#define CHARGE_EXTRA_H

#include "module_cell/unitcell.h"
#include "module_hamilt_pw/hamilt_pwdft/structure_factor.h"
#include "charge.h"

//--------------------------------
// charge extrapolation method
// written by Xiaohui 2014
//--------------------------------
class Charge_Extra
{
    public:

    Charge_Extra();
    ~Charge_Extra();

    //When Esolver is initialized, ucell.natom is not determined
    //As a result, data structures in Charge_Extra cannot be allocated
    //This is a temporary solution by delaying the allocation
    //But after ucell and Esolver are fully decoupled
    //Init_CE will be removed and everything put back in the constructor
    void Init_CE(const int& nrxx);
    void extrapolate_charge(Charge* chr, Structure_Factor* sf);
    void update_all_dis(const UnitCell& ucell);

    private:
    int istep = 0;
    int natom;
    int pot_order;
    int rho_extr;
    double omega_old;    // the old volume of the last step

    // for the second-order extrapolation
    ModuleBase::Vector3<double> *dis_old1 = nullptr;    // dis_old2 = pos_old1 - pos_old2
    ModuleBase::Vector3<double> *dis_old2 = nullptr;    // dis_old1 = pos_now - pos_old1
    ModuleBase::Vector3<double> *dis_now = nullptr;     // dis_now = pos_next - pos_now

    double** delta_rho1 = nullptr;
    double** delta_rho2 = nullptr;

    double alpha,beta;

    void find_alpha_and_beta(void);

};

#endif
