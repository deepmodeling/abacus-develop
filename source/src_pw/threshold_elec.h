#ifndef THRESHOLD_ELEC_H
#define THRESHOLD_ELEC_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_esolver/esolver_ks_lcao.h"

class Threshold_Elec
{

    friend class ModuleESolver::ESolver_KS_LCAO;
public:

    // constructor and deconstructor
    Threshold_Elec();
    ~Threshold_Elec() {};

	protected:

    double scf_thr;

    bool conv_elec;

    void set_pw_diag_thr() const;

    void update_pw_diag_thr(const int &iter);

	// this should be moved to other places, mohan note 2021-03-03
    static void print_eigenvalue(std::ofstream &ofs);
};

#endif
