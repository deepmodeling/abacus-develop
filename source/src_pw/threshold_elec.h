#ifndef THRESHOLD_ELEC_H
#define THRESHOLD_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"

class Threshold_Elec
{

	public:

    // constructor and deconstructor
    Threshold_Elec();
    ~Threshold_Elec() {};

	protected:

    double scf_thr_rho;

    bool conv_elec;

    void set_diag_thr_e() const;

    void update_diag_thr_e(const int &iter);

	// this should be moved to other places, mohan note 2021-03-03
    void print_eigenvalue(std::ofstream &ofs);
};

#endif
