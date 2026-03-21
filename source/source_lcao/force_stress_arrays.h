#ifndef FORCESTRESS_ARRAYS_H 
#define FORCESTRESS_ARRAYS_H

class ForceStressArrays
{
    public:

    ForceStressArrays(){};
    ~ForceStressArrays(){};

    //-----------------------------------------
    // force in LCAO
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DSloc_x = nullptr;
    double* DSloc_y = nullptr;
    double* DSloc_z = nullptr;

    //-----------------------------------------
    // force in LCAO
    // used in k-points algorithm.
    //-----------------------------------------
    double* DSloc_Rx = nullptr;
    double* DSloc_Ry = nullptr;
    double* DSloc_Rz = nullptr;

    //-----------------------------------------
    // dT + part of dVNL
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DHloc_fixed_x = nullptr;
    double* DHloc_fixed_y = nullptr;
    double* DHloc_fixed_z = nullptr;

    //-----------------------------------------
    // dT + part of dVNL
    // used in kpoint algorithm.
    //-----------------------------------------
    double* DHloc_fixedR_x = nullptr;
    double* DHloc_fixedR_y = nullptr;
    double* DHloc_fixedR_z = nullptr;

    //----------------------------------------
    // r_mu - r_nu
    //----------------------------------------

    double* DH_r = nullptr;//zhengdy added 2017-07

    double* stvnl11 = nullptr;
    double* stvnl12 = nullptr;
    double* stvnl13 = nullptr;
    double* stvnl22 = nullptr;
    double* stvnl23 = nullptr;
    double* stvnl33 = nullptr;

    double* DSloc_11 = nullptr;
    double* DSloc_12 = nullptr;
    double* DSloc_13 = nullptr;
    double* DSloc_22 = nullptr;
    double* DSloc_23 = nullptr;
    double* DSloc_33 = nullptr;

    double* DHloc_fixed_11 = nullptr;
    double* DHloc_fixed_12 = nullptr;
    double* DHloc_fixed_13 = nullptr;
    double* DHloc_fixed_22 = nullptr;
    double* DHloc_fixed_23 = nullptr;
    double* DHloc_fixed_33 = nullptr;

};

#endif
