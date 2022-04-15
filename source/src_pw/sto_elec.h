#ifndef STO_ELEC_H
#define STO_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "threshold_elec.h"
#include "sto_wf.h"
#include "sto_iter.h"

//----------------------------------------------
// methods based on stochastic wave functions
//----------------------------------------------

class Stochastic_Elec: private Threshold_Elec
{

	public:

    // constructor and deconstructor
    Stochastic_Elec();
    ~Stochastic_Elec();

    int iter;

    static double avg_iter;

    void scf_stochastic(const int &istep);

	private:

    Stochastic_Iter stoiter;

    void c_bands(const int &istep);

    bool check_stop_now(void);
};

#endif// Eelectrons_Stochastic
