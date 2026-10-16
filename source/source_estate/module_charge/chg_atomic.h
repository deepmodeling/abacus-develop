#ifndef CHG_ATOMIC_H
#define CHG_ATOMIC_H

#include "source_base/complexmatrix.h"
#include "source_basis/module_pw/pw_basis.h"

#include <ostream>

class UnitCell;

namespace module_charge
{

/// Configuration for atomic_rho, replacing GlobalV/PARAM reads
struct AtomicRhoCfg
{
    double nelec;        ///< target total electron number (PARAM.inp.nelec)
    int test_charge;    ///< verbosity flag (PARAM.inp.test_charge)
    bool domag;         ///< whether to compute magnetization (PARAM.globalv.domag)
    bool domag_z;       ///< whether to compute z-only magnetization
    std::ostream& ofs_warning; ///< warning output stream
};

// Superposition of atomic charges contained in the array rho_at
// (read from pseudopotential files).
void atomic_rho(const int spin_number_need,
                const double& omega,
                double** rho_in,
                const ModuleBase::ComplexMatrix& strucFac,
                const UnitCell& ucell,
                const ModulePW::PW_Basis* rhopw,
                const AtomicRhoCfg& cfg);

} // namespace module_charge

#endif // CHG_ATOMIC_H
