#ifndef ELECSATE_PRINT_H
#define ELECSATE_PRINT_H

#include "source_estate/elecstate.h"

namespace elecstate
{
    void print_format(const std::string& name, 
                    const double& value);
    
    /// @param inp the INPUT parameters whose flags decide which energy terms and
    ///        headers are printed
    /// @param two_fermi whether the run keeps two Fermi levels; derived, so it
    ///        does not live in Input_para
    void print_etot(const Magnetism& magnet,
                    const ElecState& elec,
                    const bool converged,
                    const int& iter_in,
                    const double& scf_thr,
                    const double& scf_thr_kin,
                    const double& duration,
                    const Input_para& inp,
                    const bool two_fermi,
                    const double& pw_diag_thr = 0,
                    const double& avg_iter = 0,
                    bool print = true,
                    const double& ds_rms = -1.0);
}
#endif
