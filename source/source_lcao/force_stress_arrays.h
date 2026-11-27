#ifndef FORCESTRESS_ARRAYS_H
#define FORCESTRESS_ARRAYS_H

#include <vector>

class ForceStressArrays
{
    public:

    ForceStressArrays(){};
    ~ForceStressArrays(){};

    //-----------------------------------------
    // force in LCAO
    // used in gamma only algorithm.
    //-----------------------------------------
    std::vector<double> DSloc_x;
    std::vector<double> DSloc_y;
    std::vector<double> DSloc_z;

    //-----------------------------------------
    // force in LCAO
    // used in k-points algorithm.
    //-----------------------------------------
    std::vector<double> DSloc_Rx;
    std::vector<double> DSloc_Ry;
    std::vector<double> DSloc_Rz;

    //-----------------------------------------
    // dT + part of dVNL
    // used in gamma only algorithm.
    //-----------------------------------------
    std::vector<double> DHloc_fixed_x;
    std::vector<double> DHloc_fixed_y;
    std::vector<double> DHloc_fixed_z;

    //-----------------------------------------
    // dT + part of dVNL
    // used in kpoint algorithm.
    //-----------------------------------------
    std::vector<double> DHloc_fixedR_x;
    std::vector<double> DHloc_fixedR_y;
    std::vector<double> DHloc_fixedR_z;

    //----------------------------------------
    // r_mu - r_nu
    //----------------------------------------

    std::vector<double> DH_r;//zhengdy added 2017-07

    std::vector<double> stvnl11;
    std::vector<double> stvnl12;
    std::vector<double> stvnl13;
    std::vector<double> stvnl22;
    std::vector<double> stvnl23;
    std::vector<double> stvnl33;

};

#endif
