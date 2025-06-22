#ifndef ELECSATE_PRINT_H
#define ELECSATE_PRINT_H

#include "module_elecstate/elecstate.h"

namespace elecstate
{
    void print_format(const std::string& name, 
                    const double& value);
    
    void print_etot(const Magnetism& magnet,
                    const ElecState& elec,
                    const bool converged,
                    const int& iter_in,
                    const double& scf_thr,
                    const double& scf_thr_kin,
                    const double& duration,
                    const double& pw_diag_thr = 0,
                    const double& avg_iter = 0,
                    bool print = true);
}
#endif
