#ifndef CHARGE_EXTRA_H
#define CHARGE_EXTRA_H

#include "../module_cell/unitcell_pseudo.h"

using namespace std;

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
    void Init_CE();
    void extrapolate_charge(void);
    void save_pos_next(const UnitCell_pseudo& ucell);
    void update_istep(const int &step);
    void update_all_pos(const UnitCell_pseudo& ucell);

    private:
    int istep = 0;
    int natom;
    int pot_order;
    int rho_extr;

    // for the second-order extrapolation
    ModuleBase::Vector3<double> *pos_old1;
    ModuleBase::Vector3<double> *pos_old2;
    ModuleBase::Vector3<double> *pos_now;
    ModuleBase::Vector3<double> *pos_next;

    double** delta_rho1;
    double** delta_rho2;

    double alpha,beta;

    void find_alpha_and_beta(void);

};

#endif
