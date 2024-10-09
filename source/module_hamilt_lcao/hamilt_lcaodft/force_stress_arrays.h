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
    double* DSloc_x;
    double* DSloc_y;
    double* DSloc_z;

    //-----------------------------------------
    // force in LCAO
    // used in k-points algorithm.
    //-----------------------------------------
    double* DSloc_Rx;
    double* DSloc_Ry;
    double* DSloc_Rz;

    //-----------------------------------------
    // dT + part of dVNL
    // used in gamma only algorithm.
    //-----------------------------------------
    double* DHloc_fixed_x;
    double* DHloc_fixed_y;
    double* DHloc_fixed_z;

    //-----------------------------------------
    // dT + part of dVNL
    // used in kpoint algorithm.
    //-----------------------------------------
    double* DHloc_fixedR_x;
    double* DHloc_fixedR_y;
    double* DHloc_fixedR_z;

    //----------------------------------------
    // r_mu - r_nu
    //----------------------------------------

    double* DH_r;//zhengdy added 2017-07

    double* stvnl11;
    double* stvnl12;
    double* stvnl13;
    double* stvnl22;
    double* stvnl23;
    double* stvnl33;

    double* DSloc_11;
    double* DSloc_12;
    double* DSloc_13;
    double* DSloc_22;
    double* DSloc_23;
    double* DSloc_33;

    double* DHloc_fixed_11;
    double* DHloc_fixed_12;
    double* DHloc_fixed_13;
    double* DHloc_fixed_22;
    double* DHloc_fixed_23;
    double* DHloc_fixed_33;

};

#endif
