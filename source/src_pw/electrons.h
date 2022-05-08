#ifndef INCLUDE_ELECTRONS_H
#define INCLUDE_ELECTRONS_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "threshold_elec.h"

class Electrons: private Threshold_Elec
{

public:

    // constructor and deconstructor
    Electrons();
    ~Electrons();

public:

    //==========================================================
    // avg_iter : Average iteration number in ccgdiagg
    //==========================================================

    int iter;

    static double avg_iter;

    int test;

    int unit;

    double delta_total_energy;

    // It should be replaced by ESolver_KS::Run. GPU still use it temperarily.
    void self_consistent(const int &istep);

    void non_self_consistent(const int &istep);

    void c_bands(const int &istep);

private:
    bool check_stop_now(void);

    void init_mixstep_final_scf(void);
};

#endif
